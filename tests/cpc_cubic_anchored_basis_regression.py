#!/usr/bin/env python3
"""Prove exact ED2 cubic translation-anchored basis states and normalizer bits.

Execute generic atomic-orbits vs anchored-cubic basis on identical input decks,
compare every sorted representative rank and IEEE-754 list_r bit, and check
64/216-site full native Hamiltonian/Lanczos outputs. Python >=3.6, stdlib.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from cpc_cubic_pointgroup import canonical_bonds, sector_dimensions, symmetry_generators, write_input
from cpc_pointgroup_profile import run_case


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while True:
            b = stream.read(1048576)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--driver', required=True, type=Path)
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--sides', default='4,6,8,10')
    p.add_argument('--threads', default='1,4')
    p.add_argument('--sectors', default='oh-a1g,oh-a2g,oh-a1u,oh-a2u')
    p.add_argument('--timeout', type=int, default=600)
    p.add_argument('--skip-lanczos', action='store_true')
    args = p.parse_args()
    driver = args.driver.resolve(strict=True)
    exe = args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix='ed2-cubic-anchored-basis-'))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    sectors = args.sectors.split(',')
    assert sectors and set(sectors) <= {'oh-a1g', 'oh-a2g', 'oh-a1u', 'oh-a2u'}
    report = []
    for side in (int(v) for v in args.sides.split(',')):
        dimensions = sector_dimensions(side)
        bonds = canonical_bonds(side)
        for sector in sectors:
            for threads in (int(v) for v in args.threads.split(',')):
                captures = []
                for anchored in (False, True):
                    work = root / 's{}-{}-t{}-{}'.format(
                        side, sector, threads, 'anchored' if anchored else 'old')
                    write_input(work, side, 3, symmetry_generators(side, sector), bonds, mnte=64)
                    env = dict(os.environ)
                    env.update(ED2_CUBIC_FAST_GROUP_VALIDATION='1',
                               ED2_CUBIC_BASIS_PROOF_DUMP='1',
                               ED2_BASIS_DETAIL_TIMING='1',
                               ED2_ALLOW_LARGE_GROUP_VALIDATION='1',
                               ED2_CUBIC_ANCHORED_REP='1',
                               OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
                               OPENBLAS_NUM_THREADS='1')
                    env.pop('ED2_FORCE_LEGACY_ORBITS', None)
                    if anchored:
                        env['ED2_CUBIC_ANCHORED_BASIS'] = '1'
                    else:
                        env.pop('ED2_CUBIC_ANCHORED_BASIS', None)
                    with (work / 'input.dat').open('rb') as inp, (work / 'driver_stdout.log').open('wb') as out:
                        result = subprocess.run([str(driver)], stdin=inp, stdout=out,
                                                stderr=subprocess.STDOUT, cwd=work,
                                                env=env, timeout=args.timeout, check=False)
                    log = (work / 'driver_stdout.log').read_text(errors='replace')
                    assert result.returncode == 0, (work, result.returncode, log[-2500:])
                    assert ('Basis builder: anchored-cubic' in log) is anchored
                    assert 'Symmetry certificate: canonical cubic O_h semidirect product verified' in log
                    d = re.findall(r'ED2_BASIS_DRIVER_DIM=(\d+)', log)
                    assert len(d) == 1 and int(d[0]) == dimensions[sector]
                    bits = (work / 'cubic_basis_proof_dump.txt').read_bytes()
                    assert len(bits.splitlines()) == dimensions[sector]
                    times = re.findall(r'ED2_BASIS_PHASE orbit_scan wall_s=\s*([+\-0-9.eE]+)', log)
                    assert len(times) == 1 and float(times[0]) >= 0
                    captures.append((digest(work / 'input.dat'), bits, float(times[0])))
                assert captures[0][0] == captures[1][0]
                assert captures[0][1] == captures[1][1], ('EXACT BITWISE BASES MISMATCH', side, sector, threads)
                report.append(dict(side=side,sector=sector,threads=threads,
                                   dimension=dimensions[sector],
                                   exact_rank_and_norm_bits=True,
                                   original_orbit_s=captures[0][2],
                                   anchored_orbit_s=captures[1][2],
                                   proof_sha256=hashlib.sha256(captures[1][1]).hexdigest()))
                print('ANCHORED_BASIS_BITWISE_PASS side={} sector={} threads={} D={} '
                      'original={:.6f}s new={:.6f}s'.format(
                          side, sector, threads, dimensions[sector],
                          captures[0][2], captures[1][2]),flush=True)
                if side in (4,6) and threads == 1 and not args.skip_lanczos:
                    physical = []
                    for anchored in (False, True):
                        work = root / 'physical-s{}-{}-{}'.format(
                            side, sector, 'anchored' if anchored else 'original')
                        write_input(work, side, 3, symmetry_generators(side, sector), bonds, mnte=64)
                        keys=('ED2_CUBIC_ANCHORED_BASIS','ED2_CUBIC_FAST_GROUP_VALIDATION',
                              'ED2_CUBIC_ANCHORED_REP')
                        prior = {name:os.environ.get(name) for name in keys}
                        try:
                            os.environ['ED2_CUBIC_FAST_GROUP_VALIDATION']='1'
                            os.environ['ED2_CUBIC_ANCHORED_REP']='1'
                            if anchored:
                                os.environ['ED2_CUBIC_ANCHORED_BASIS']='1'
                            else:
                                os.environ.pop('ED2_CUBIC_ANCHORED_BASIS',None)
                            physical.append(run_case(exe,work,side**3,3,sector,
                                                     dimensions[sector],1,'two-stage',args.timeout))
                        finally:
                            for name,value in prior.items():
                                if value is None:
                                    os.environ.pop(name,None)
                                else:
                                    os.environ[name]=value
                    assert physical[0]['input_sha256']==physical[1]['input_sha256']
                    assert abs(physical[0]['energy']-physical[1]['energy'])<1e-8
                    assert physical[0]['residual']<=1.01e-8 and physical[1]['residual']<=1.01e-8
                    assert physical[0]['h_applications']==physical[1]['h_applications']
                    print('ANCHORED_BASIS_FULL_HAMILTONIAN_PASS side={} sector={}'.format(
                        side,sector),flush=True)
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('CUBIC_ANCHORED_BASIS_REGRESSION_PASS {} exact bitwise basis pairs'.format(len(report)),
          flush=True)


if __name__=='__main__':
    main()
