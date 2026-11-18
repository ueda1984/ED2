#!/usr/bin/env python3
"""1000-site full cubic O_h normal-form smoke at fixed one-down occupancy.

This verifies large-group validation and a trivially solvable one-magnon
physical energy. It must NEVER be reported as a three-magnon benchmark.
"""
import argparse
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

from cpc_cubic_pointgroup import canonical_bonds, symmetry_generators, write_input


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    exe=a.exe.resolve(strict=True)
    root=a.output if a.output else Path(tempfile.mkdtemp(prefix='ed2-cubic-oh-1000-n1-'))
    side=10
    work=root if a.output else root/'single'
    write_input(work,side,1,symmetry_generators(side,'oh-a1g'),
                canonical_bonds(side),400,-1)
    env=dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
             OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
             ED2_LANCZOS_TWO_STAGE='1', ED2_ALLOW_LARGE_GROUP_VALIDATION='1',
             ED2_PHASE_TIMING='1', ED2_DIAGNOSTICS='diagnostics.txt', ED2_RNG_SEED='1984')
    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as output:
        proc=subprocess.run([str(exe)],cwd=work,stdin=inp,stdout=output,
                            stderr=subprocess.STDOUT,env=env,timeout=180)
    log=(work/'stdout.log').read_text(errors='replace')
    assert proc.returncode==0, log[-2500:]
    assert 'Symmetry action: noncommuting normal form and character validated' in log
    assert 'Basis builder: atomic-orbits' in log
    assert re.search(r'THS\(k\)\s*=\s*1\b',log)
    energy=float((work/'eigenvalues.dat').read_text().splitlines()[0].split()[-1])
    # Uniform one-down state: vacuum ZZ + onsite ZZ shift + XY hopping.
    expected=(3*side**3*.77+6*side**3*.19)/4 - (6*.77+12*.19)/2 + (6+12*.31)/2
    assert math.isclose(expected, 863.91,abs_tol=1e-10)
    assert abs(energy-expected) < 1e-8, (energy,expected)
    lines=(work/'diagnostics.txt').read_text().splitlines()
    rows=[v for s in lines if (v:=s.split()) and len(v)==8 and v[0]=='1' and v[1:3]==['T','T']]
    assert len(rows)==1
    assert float(rows[0][5]) <= 1e-8
    print(f'PASS cubic O_h 1000-site one-magnon fast path: D=1 E={energy:.12f} residual={float(rows[0][5]):.3e}',flush=True)

if __name__=='__main__':
    main()
