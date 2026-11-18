#!/usr/bin/env python3
"""Compare every native cubic cache hop against legacy symmetry ID and phase.

Only available for the exact Gamma/O_h, spin-1/2, 3-down J1-J2 XXZ cube.
Both the whole spectrum reference (small cube) and fixed-seed Lanczos
checks are independent of the new lookup. No performance acceptance gate.
"""
import argparse
import json
import os
from pathlib import Path
import tempfile

from cpc_cubic_pointgroup import (SECTORS, canonical_bonds, sector_dimensions,
                                  symmetry_generators, write_input)
from cpc_pointgroup_profile import run_case


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--sizes', default='4,6')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--verify-hops', action='store_true')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix='ed2-native-cubic-anchor-'))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    exe = args.exe.resolve(strict=True)
    results = []
    for side in (int(s) for s in args.sizes.split(',')):
        dimensions = sector_dimensions(side)
        bonds = canonical_bonds(side)
        for sector in SECTORS[1:]:
            reference = {}
            for mode in ('legacy', 'anchored'):
                work = root / ('s{}-{}-{}-t{}'.format(side,sector,mode,args.threads))
                write_input(work, side, 3, symmetry_generators(side,sector), bonds,
                            mnte=64)
                old_anchor = os.environ.get('ED2_CUBIC_ANCHORED_REP')
                old_verify = os.environ.get('ED2_CUBIC_ANCHORED_VERIFY')
                try:
                    if mode == 'anchored':
                        os.environ['ED2_CUBIC_ANCHORED_REP'] = '1'
                        if args.verify_hops:
                            os.environ['ED2_CUBIC_ANCHORED_VERIFY'] = '1'
                        else:
                            os.environ.pop('ED2_CUBIC_ANCHORED_VERIFY', None)
                    else:
                        os.environ.pop('ED2_CUBIC_ANCHORED_REP', None)
                        os.environ.pop('ED2_CUBIC_ANCHORED_VERIFY', None)
                    obs = run_case(exe,work,side**3,3,sector,dimensions[sector],
                                   args.threads,'two-stage',600)
                finally:
                    for name, value in [('ED2_CUBIC_ANCHORED_REP',old_anchor),
                                        ('ED2_CUBIC_ANCHORED_VERIFY',old_verify)]:
                        if value is None:
                            os.environ.pop(name,None)
                        else:
                            os.environ[name] = value
                stdout=(work/'stdout.log').read_text(errors='replace')
                assert ('Representative lookup: cubic-anchored-orbits' in stdout) == (mode=='anchored')
                if mode == 'anchored' and args.verify_hops:
                    assert 'per-hop legacy verification=T' in stdout
                reference[mode] = obs
                print('PASS s{} {} {} t{} E={:.14f} R={:.3e} cache={:.6f}s'.format(
                    side,sector,mode,args.threads,obs['energy'],obs['residual'],
                    obs['phase_wall_s']['operator_cache']),flush=True)
            a,b=reference['legacy'],reference['anchored']
            assert a['input_sha256'] == b['input_sha256']
            assert a['dimension']==b['dimension']
            assert abs(a['energy']-b['energy'])<=1e-8
            assert a['residual'] <= 1.01e-8 and b['residual']<=1.01e-8
            # At OMP>1, floating-point summation and Krylov termination can
            # change iteration counts despite identical matrix elements.
            if args.threads == 1:
                assert a['h_applications'] == b['h_applications']
            results.append(dict(side=side,sector=sector,threads=args.threads,
                                energy=a['energy'],direct_residual=b['residual'],
                                basis=a['dimension'],
                                legacy_h_applications=a['h_applications'],
                                anchored_h_applications=b['h_applications'],
                                legacy_cache_s=a['phase_wall_s']['operator_cache'],
                                anchored_cache_s=b['phase_wall_s']['operator_cache'],
                                verify_hops=args.verify_hops))
    (root/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    print('CUBIC_ANCHORED_NATIVE_REGRESSION_PASS {} matched sectors (verify_hops={})'.format(
        len(results),args.verify_hops),flush=True)


if __name__ == '__main__':
    main()
