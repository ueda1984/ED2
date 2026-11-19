#!/usr/bin/env python3
"""Check Lanczos scalar-shift rejection and ordinary resolved/degenerate cases.

The three-site one-deviation model has H = C I + H0, where H0 has
eigenvalues (-0.005, 0, 0.005). At C=1e12 the gap is about 41 binary64
ULPs. A small residual relative to C alone must not certify closure.
Uses only the Python standard library and the supplied native executable.
"""
import argparse
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys


DIAGNOSTIC = 'Lanczos numerical near-breakdown cannot certify an invariant subspace'
RESULT_FILES = ('eigenvalues.dat', 'local_mag.dat', 'two_body_cf_xyz.dat',
                'two_body_cf_z+-.dat', 'energy.dat', 'wf1.bin')


def prepare(work, name):
    work.mkdir()
    if name == 'one-dimensional':
        n, lo, hi, dim = 1, 1, 1, 1
        fields, bonds = '1 0 0 1.4\n', ''
        expected_energy, expected_mag = -.7, [-.5]
    elif name in ('centered', 'large-shift', 'scalar-degenerate'):
        n, lo, hi, dim = 3, 1, 1, 3
        shift = 1e12 if name == 'large-shift' else (1.0 if name == 'scalar-degenerate' else 0.0)
        fields = ''.join(f'{i} 0 0 {2 * shift:.17g}\n' for i in range(1, n + 1))
        bonds = '' if name == 'scalar-degenerate' else '1 2 .01 .01 0 0 0 0 0 0 0\n'
        expected_energy = shift if name == 'scalar-degenerate' else shift - .005
        expected_mag = None if name == 'scalar-degenerate' else [0., 0., .5]
    else:
        n, lo, hi, dim = 2, 0, 2, 4
        fields = ''
        bonds = '1 2 1 1 1 0 0 0 0 0 0\n' if name == 'heisenberg-dimer' else ''
        expected_energy = -.75 if name == 'heisenberg-dimer' else 0.0
        expected_mag = [0., 0.] if name == 'heisenberg-dimer' else None
    (work / 'spins.dat').write_text('0.5\n' * n)
    (work / 'caps.dat').write_text('1\n' * n)
    (work / 'fields.dat').write_text(fields)
    (work / 'bonds.dat').write_text(bonds)
    (work / 'pairs.dat').write_text('1 2\n' if n > 1 else '')
    (work / 'input.dat').write_text(f'''&input_parameters
NOS={n},NODmin={lo},NODmax={hi},ALG=1,
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={len(fields.splitlines())},NO_two={len(bonds.splitlines())},wk_dim=-1,MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=1,re_wf=0,cal_lm=1,
cal_cf={int(n > 1)},NO_two_cf={int(n > 1)},FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-14,minitr=1,maxitr=100,itrint=1
/
''')
    return dim, expected_energy, expected_mag


def norm_of_native_checkpoint(path, dim):
    data = path.read_bytes()
    for marker in (4, 8):
        if (len(data) == dim * 16 + 2 * marker
                and int.from_bytes(data[:marker], sys.byteorder, signed=True) == dim * 16
                and data[:marker] == data[-marker:]):
            numbers = struct.unpack(f'={2 * dim}d', data[marker:-marker])
            assert all(math.isfinite(x) for x in numbers)
            return sum(x * x for x in numbers)
    raise AssertionError('Unrecognized native checkpoint record layout')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    exe, root = args.exe.resolve(), args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = {k: v for k, v in os.environ.items() if not k.startswith('ED2_')}
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    results = []
    for name in ('zero-hamiltonian', 'one-dimensional', 'centered',
                 'heisenberg-dimer', 'large-shift', 'scalar-degenerate'):
        work = root / name
        dim, expected_energy, expected_mag = prepare(work, name)
        before = {}
        if name == 'large-shift':
            for filename in RESULT_FILES:
                (work / filename).write_bytes(b'previous completed result\n')
            before = {p.name: p.read_bytes() for p in work.iterdir() if p.is_file()}
        with (work / 'input.dat').open() as inp, (work / 'stdout.log').open('w') as out:
            proc = subprocess.run([str(exe)], cwd=work, env=env, stdin=inp,
                                  stdout=out, stderr=subprocess.STDOUT, timeout=30)
        log = (work / 'stdout.log').read_text(errors='replace')
        result = dict(name=name, returncode=proc.returncode, threads=args.threads)
        try:
            if name == 'large-shift':
                assert proc.returncode != 0, 'A numerically unresolved shifted model was accepted'
                assert DIAGNOSTIC in log, 'Missing numerical near-breakdown diagnostic'
                assert 'one-dimensional invariant Krylov subspace' not in log
                for filename, data in before.items():
                    assert (work / filename).read_bytes() == data, f'Failure modified {filename}'
                extra = {p.name for p in work.iterdir()} - before.keys() - {'stdout.log'}
                assert not extra, f'Failure created files: {sorted(extra)}'
                result['behavior'] = 'controlled rejection; prior outputs and checkpoint unchanged'
            elif name == 'scalar-degenerate' and proc.returncode != 0:
                # With H=I a floating-point residual may be exactly zero or
                # a few ULPs. Without an observed spectral width the latter
                # is deliberately rejected, even though this model is valid.
                assert DIAGNOSTIC in log
                assert not any((work / file).exists() for file in RESULT_FILES)
                result['behavior'] = 'documented conservative rejection of scalar degeneracy'
            else:
                assert proc.returncode == 0, 'Ordinary resolved/degenerate model was rejected'
                energy = float((work / 'eigenvalues.dat').read_text().split()[1])
                assert math.isfinite(energy)
                assert abs(energy - expected_energy) < 2e-12, (energy, expected_energy)
                norm2 = norm_of_native_checkpoint(work / 'wf1.bin', dim)
                assert abs(norm2 - 1) < 2e-12, norm2
                rows = [[float(v) for v in line.split()] for line in
                        (work / 'local_mag.dat').read_text().splitlines()]
                assert all(math.isfinite(v) for row in rows for v in row)
                assert all(abs(v) <= .5 + 2e-12 for row in rows for v in row[1:])
                if expected_mag is not None:
                    assert max(abs(row[3] - value) for row, value in zip(rows, expected_mag)) < 2e-11
                result.update(energy=energy, norm2=norm2, behavior='finite normalized analytic eigenstate')
            result['status'] = 'PASS'
        except Exception as error:
            result.update(status='FAIL', error=str(error))
        results.append(result)
    summary = dict(executable=str(exe), results=results,
                   passed=sum(row['status'] == 'PASS' for row in results),
                   failed=sum(row['status'] == 'FAIL' for row in results))
    (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return int(summary['failed'] != 0)


if __name__ == '__main__':
    raise SystemExit(main())
