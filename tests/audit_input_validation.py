#!/usr/bin/env python3
"""Native input-range, finite-input and read-only restart rejection regression.

Run against both an -O2 build and GNU's default -Ofast build: floating-point
optimizations must not remove the NaN/Inf checks. Uses the standard library.
Checkpoints are produced by the executable itself, retaining its record format.
"""
import argparse
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys


OUTPUTS = ('eigenvalues.dat', 'local_mag.dat', 'two_body_cf_xyz.dat',
           'two_body_cf_z+-.dat')


def prepare(work, *, algorithm=3, restart=False, tolerance='1d-12'):
    work.mkdir()
    (work / 'spins.dat').write_text('0.5\n0.5\n')
    (work / 'caps.dat').write_text('1\n1\n')
    (work / 'fields.dat').write_text('1 0.2 0.3 1\n2 -0.1 0.4 0.7\n')
    (work / 'bonds.dat').write_text('1 2 0.3 0.3 0.3 0 0 0 0 0 0\n')
    (work / 'pairs.dat').write_text('1 2\n')
    (work / 'input.dat').write_text(f'''&input_parameters
NOS=2,NODmin=0,NODmax=2,ALG={algorithm},
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one=2,NO_two=1,wk_dim=-1,MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=1,re_wf={int(restart)},
cal_lm=1,cal_cf=1,NO_two_cf=1,FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv={tolerance},minitr=1,maxitr=100,itrint=1
/
&input_TRLan
NOE=1,NOK=1,NOM=3,maxitr=500,lnc_ene_conv={tolerance},i_vec_min=1,i_vec_max=1
/
''')


def run(exe, work):
    env = {k: v for k, v in os.environ.items() if not k.startswith('ED2_')}
    env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (work / 'input.dat').open() as inp, (work / 'stdout.log').open('w') as out:
        proc = subprocess.run([str(exe)], cwd=work, env=env, stdin=inp,
                              stdout=out, stderr=subprocess.STDOUT, timeout=30)
    return proc.returncode, (work / 'stdout.log').read_text(errors='replace')


def read_numbers(path):
    return [float(x.replace('D', 'E')) for x in path.read_text().split()]


def checkpoint_marker_size(data):
    # GNU and Intel default native sequential unformatted records use a
    # four-byte marker; accept eight-byte markers without fabricating a file.
    payload_size = 4 * 16  # Four complex binary64 coefficients in this model.
    for size in (4, 8):
        if (len(data) == payload_size + 2 * size
                and int.from_bytes(data[:size], sys.byteorder, signed=True) == payload_size
                and data[-size:] == data[:size]):
            return size
    raise AssertionError('Unrecognized native checkpoint record layout')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    exe, root = args.exe.resolve(), args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    results = []
    seed = root / 'valid-finite-dense'
    prepare(seed)
    status, _ = run(exe, seed)
    if status != 0:
        raise RuntimeError(f'Finite control failed; see {seed}/stdout.log')
    for filename in OUTPUTS:
        values = read_numbers(seed / filename)
        if not values or not all(math.isfinite(x) for x in values):
            raise AssertionError(f'Finite control produced invalid {filename}')
    wave = (seed / 'wf1.bin').read_bytes()
    marker = checkpoint_marker_size(wave)
    energy = (seed / 'energy.dat').read_bytes()
    results.append(dict(name=seed.name, status='PASS', returncode=status))

    def checkpoint(work):
        (work / 'wf1.bin').write_bytes(wave)
        (work / 'energy.dat').write_bytes(energy)

    def reject(work, diagnostic):
        # A failed invocation must also leave previous scientific outputs
        # untouched, including when wr_wf=1 is present in a restart deck.
        for filename in OUTPUTS:
            (work / filename).write_bytes(b'previous completed result\n')
        before = {p.name: p.read_bytes() for p in work.iterdir() if p.is_file()}
        status, log = run(exe, work)
        error = None
        if status == 0:
            error = 'Invalid input returned success'
        elif diagnostic and diagnostic not in log:
            error = f'Missing diagnostic: {diagnostic}'
        elif not diagnostic and not any(word in log.lower() for word in ('error', 'severe')):
            error = 'Malformed restart produced no clear error diagnostic'
        for filename, content in before.items():
            if not (work / filename).is_file() or (work / filename).read_bytes() != content:
                error = f'Rejected input changed {filename}'
        after_names = {p.name for p in work.iterdir() if p.is_file()}
        unexpected = after_names - before.keys() - {'stdout.log'}
        if unexpected:
            error = f'Rejected input created unexpected files: {sorted(unexpected)}'
        result = dict(name=work.name, returncode=status, status='FAIL' if error else 'PASS')
        if error:
            result['error'] = error
        results.append(result)

    # These read the genuine normalized checkpoint in both solver branches.
    for algorithm in (1, 2):
        work = root / f'valid-finite-restart-alg{algorithm}'
        prepare(work, algorithm=algorithm, restart=True)
        checkpoint(work)
        status, _ = run(exe, work)
        error = None
        if status != 0:
            error = 'Finite checkpoint was rejected'
        elif (work / 'wf1.bin').read_bytes() != wave or (work / 'energy.dat').read_bytes() != energy:
            error = 'Restart modified its checkpoint'
        elif not all(math.isfinite(x) for filename in OUTPUTS[1:]
                     for x in read_numbers(work / filename)):
            error = 'Finite restart produced nonfinite observables'
        results.append(dict(name=work.name, returncode=status,
                            status='FAIL' if error else 'PASS', error=error))

    # Include negative zero and subnormal input in the accepted finite domain.
    work = root / 'valid-finite-zero-and-subnormal'
    prepare(work)
    (work / 'fields.dat').write_text('1 -0.0 4.9406564584124654e-324 1\n2 -0.1 0.4 0.7\n')
    status, _ = run(exe, work)
    valid = status == 0 and all(math.isfinite(x) for x in read_numbers(work / 'eigenvalues.dat'))
    results.append(dict(name=work.name, returncode=status, status='PASS' if valid else 'FAIL'))

    # Valid edge values stay accepted: S=0/cap=0, integer spin/full cap,
    # and the largest half-integer safely convertible to default INTEGER.
    for label, spins, caps, max_nod, no_interactions in (
        ('zero-spin', '0\n0.5\n', '0\n1\n', 1, True),
        ('integer-spin', '1\n0.5\n', '2\n1\n', 3, False),
        ('spin-conversion-bound', '1073741823.5\n0\n', '0\n0\n', 0, True),
    ):
        work = root / f'valid-{label}'
        prepare(work)
        (work / 'spins.dat').write_text(spins)
        (work / 'caps.dat').write_text(caps)
        deck = (work / 'input.dat').read_text().replace('NODmax=2', f'NODmax={max_nod}')
        if no_interactions:
            deck = deck.replace('NO_one=2,NO_two=1', 'NO_one=0,NO_two=0')
            (work / 'fields.dat').write_text('')
            (work / 'bonds.dat').write_text('')
        (work / 'input.dat').write_text(deck)
        status, _ = run(exe, work)
        valid = status == 0 and all(math.isfinite(x) for x in read_numbers(work / 'eigenvalues.dat'))
        results.append(dict(name=work.name, returncode=status, status='PASS' if valid else 'FAIL'))

    for label, old, new, diagnostic in (
        ('algorithm-zero', 'ALG=3', 'ALG=0', 'ALG must be 1, 2, or 3'),
        ('algorithm-large', 'ALG=3', 'ALG=4', 'ALG must be 1, 2, or 3'),
        ('sites-zero', 'NOS=2', 'NOS=0', 'NOS must be positive'),
        ('sites-negative', 'NOS=2', 'NOS=-1', 'NOS must be positive'),
        ('one-body-count-negative', 'NO_one=2', 'NO_one=-1', 'Interaction counts must be nonnegative'),
        ('two-body-count-negative', 'NO_two=1', 'NO_two=-1', 'Interaction counts must be nonnegative'),
        ('period-zero', 'L1=1', 'L1=0', 'Symmetry periods must be positive'),
        ('period-negative', 'L6=1', 'L6=-1', 'Symmetry periods must be positive'),
        ('nod-negative', 'NODmin=0', 'NODmin=-1', 'Require 0 <= NODmin <= NODmax'),
        ('nod-reversed', 'NODmin=0', 'NODmin=3', 'Require 0 <= NODmin <= NODmax'),
        ('correlation-count-negative', 'NO_two_cf=1', 'NO_two_cf=-1', 'NO_two_cf must be nonnegative'),
    ):
        work = root / label
        prepare(work)
        deck = (work / 'input.dat').read_text()
        assert deck.count(old) == 1
        (work / 'input.dat').write_text(deck.replace(old, new))
        reject(work, diagnostic)

    for label, filename, contents, diagnostic in (
        ('field-site-zero', 'fields.dat', '0 0 0 1\n2 0 0 1\n', 'Field index must be a site in 1..NOS'),
        ('field-site-large', 'fields.dat', '3 0 0 1\n2 0 0 1\n', 'Field index must be a site in 1..NOS'),
        ('bond-site-zero', 'bonds.dat', '0 2 1 1 1 0 0 0 0 0 0\n', 'Bond indices must be distinct sites in 1..NOS'),
        ('bond-site-large', 'bonds.dat', '1 3 1 1 1 0 0 0 0 0 0\n', 'Bond indices must be distinct sites in 1..NOS'),
        ('bond-self', 'bonds.dat', '1 1 1 1 1 0 0 0 0 0 0\n', 'Bond indices must be distinct sites in 1..NOS'),
        ('correlation-site-zero', 'pairs.dat', '0 2\n', 'Correlation indices must be distinct sites in 1..NOS'),
        ('correlation-site-large', 'pairs.dat', '1 3\n', 'Correlation indices must be distinct sites in 1..NOS'),
        ('correlation-self', 'pairs.dat', '1 1\n', 'Correlation indices must be distinct sites in 1..NOS'),
        ('spin-negative', 'spins.dat', '-0.5\n0.5\n', 'Spin must be nonnegative'),
        ('spin-fractional', 'spins.dat', '0.7\n0.5\n', 'Spin must be an integer or half-integer'),
        ('spin-conversion-overflow', 'spins.dat', '1073741824\n0.5\n', '2*S within the default INTEGER range'),
        ('spin-huge-finite', 'spins.dat', '1e300\n0.5\n', '2*S within the default INTEGER range'),
        ('cap-negative', 'caps.dat', '-1\n1\n', 'Local occupation cap must satisfy 0 <= cap <= 2*S'),
        ('cap-exceeds-spin', 'caps.dat', '2\n1\n', 'Local occupation cap must satisfy 0 <= cap <= 2*S'),
    ):
        work = root / label
        prepare(work)
        (work / filename).write_text(contents)
        reject(work, diagnostic)

    for label, token, value in (('nan', 'NaN', float('nan')),
                                ('positive-inf', 'Infinity', float('inf')),
                                ('negative-inf', '-Infinity', -float('inf'))):
        for target, diagnostic in (
            ('field', 'FILE_hvec contains a nonfinite field'),
            ('coupling', 'FILE_xyz_dm_gamma contains a nonfinite coupling'),
            ('spin', 'FILE_spin contains a nonfinite spin'),
        ):
            work = root / f'{target}-{label}'
            prepare(work)
            if target == 'field':
                (work / 'fields.dat').write_text(f'1 0 0 {token}\n2 -0.1 0.4 0.7\n')
            elif target == 'coupling':
                (work / 'bonds.dat').write_text(f'1 2 0.3 0.3 0.3 0 0 0 0 0 {token}\n')
            else:
                (work / 'spins.dat').write_text(f'0.5\n{token}\n')
            reject(work, diagnostic)
        for algorithm in (1, 2):
            for restart in (False, True):
                work = root / f'tolerance-alg{algorithm}-restart{int(restart)}-{label}'
                prepare(work, algorithm=algorithm, restart=restart, tolerance=token)
                if restart:
                    checkpoint(work)
                reject(work, 'tolerance lnc_ene_conv must be finite')
        for component in ('real', 'imaginary'):
            work = root / f'checkpoint-{component}-{label}'
            prepare(work, algorithm=1, restart=True)
            checkpoint(work)
            data = bytearray(wave)
            offset = marker + (0 if component == 'real' else 8)
            data[offset:offset + 8] = struct.pack('=d', value)
            (work / 'wf1.bin').write_bytes(data)
            reject(work, 'Saved wavefunction contains a nonfinite coefficient')
        work = root / f'checkpoint-energy-{label}'
        prepare(work, algorithm=2, restart=True)
        checkpoint(work)
        (work / 'energy.dat').write_text(f'1 {token}\n')
        reject(work, 'Saved energy.dat contains a nonfinite energy')

    for kind in ('truncated-wavefunction', 'missing-wavefunction', 'malformed-energy'):
        work = root / f'checkpoint-{kind}'
        prepare(work, algorithm=1, restart=True)
        checkpoint(work)
        if kind == 'truncated-wavefunction':
            (work / 'wf1.bin').write_bytes(wave[:marker + 8])
        elif kind == 'missing-wavefunction':
            (work / 'wf1.bin').unlink()
        else:
            (work / 'energy.dat').write_text('1 not-an-energy\n')
        reject(work, None)

    summary = dict(executable=str(exe), results=results,
                   passed=sum(row['status'] == 'PASS' for row in results),
                   failed=sum(row['status'] == 'FAIL' for row in results))
    (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return int(summary['failed'] != 0)


if __name__ == '__main__':
    raise SystemExit(main())
