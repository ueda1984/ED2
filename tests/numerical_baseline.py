#!/usr/bin/env python3
"""Black-box audit of the unchanged ED2 executable; NumPy is test-only.

Run: python3 tests/numerical_baseline.py --exe /path/to/QS3ED2 --output /outside/checkout/results
All input decks, stdout, exit statuses and numerical assertions are retained.
The command returns nonzero for any failed case, including known defects.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

ATOL = 1e-9
RTOL = 1e-9


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model(n, bonds, fields):
    """Independent product-basis spin-1/2 operator construction, up before down."""
    pauli = [np.array([[0, 1], [1, 0]]) / 2,
             np.array([[0, -1j], [1j, 0]]) / 2,
             np.diag([.5, -.5])]
    spin = []
    for site in range(n):
        axes = []
        for op in pauli:
            a = np.ones((1, 1))
            for j in range(n):
                a = np.kron(a, op if j == site else np.eye(2))
            axes.append(a)
        spin.append(axes)
    h = np.zeros((2**n, 2**n), dtype=complex)
    for i, field in enumerate(fields):
        for a in range(3):
            h += field[a] * spin[i][a]
    for i, j, c in bonds:
        for a in range(3):
            h += c[a] * spin[i][a] @ spin[j][a]
        for a, b, k in [(1, 2, 0), (2, 0, 1), (0, 1, 2)]:
            h += c[3+k] * (spin[i][a] @ spin[j][b] - spin[i][b] @ spin[j][a])
            h += c[6+k] * (spin[i][a] @ spin[j][b] + spin[i][b] @ spin[j][a])
    np.testing.assert_allclose(h, h.conj().T, atol=1e-14, rtol=0)
    return h, spin


def projection(n, lo, hi, momentum):
    states = [i for i in range(2**n) if lo <= i.bit_count() <= hi]
    if momentum is None:
        return np.eye(2**n)[:, states]
    # Sum cyclic translations with exp(-ikr); SVD removes zero/dependent orbits.
    t = np.zeros((2**n, 2**n), complex)
    for i in range(2**n):
        j = ((i << 1) & (2**n - 1)) | (i >> (n-1))
        t[j, i] = 1
    projector = sum(np.exp(-2j*np.pi*momentum*r/n) * np.linalg.matrix_power(t, r)
                    for r in range(n)) / n
    u, singular, _ = np.linalg.svd(projector[:, states], full_matrices=False)
    return u[:, singular > 1e-10]


def run_case(exe, root, name, n, bonds, fields, alg=3, lo=0, hi=None,
             momentum=None, cache=-1, threads=1, observables=True, minitr=5):
    hi = n if hi is None else hi
    work = root / name
    work.mkdir()  # Refuse to reuse stale output directories.
    h, spin = model(n, bonds, fields)
    p = projection(n, lo, hi, momentum)
    # Truncation means PHP; invariance outside that subspace is not assumed.
    values, vectors = np.linalg.eigh(p.conj().T @ h @ p)
    ground = p @ vectors[:, 0]
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    (work/'bonds.dat').write_text(''.join(f'{i+1} {j+1} ' + ' '.join(map(str,c)) + '\n'
                                           for i,j,c in bonds))
    (work/'fields.dat').write_text(''.join(f'{i+1} ' + ' '.join(map(str,f)) + '\n'
                                            for i,f in enumerate(fields)))
    (work/'spins.dat').write_text('0.5\n'*n)
    (work/'caps.dat').write_text('1\n'*n)
    (work/'pairs.dat').write_text(''.join(f'{i+1} {j+1}\n' for i,j in pairs))
    (work/'shift.dat').write_text(''.join(f'{(i+1)%n+1}\n' for i in range(n)))
    text = f'''&input_parameters
NOS={n}, NODmin={lo}, NODmax={hi}, ALG={alg},
L1={n if momentum is not None else 1}, L2=1,L3=1,L4=1,L5=1,L6=1,
M1={momentum or 0},M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={n},NO_two={len(bonds)},wk_dim={cache},MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',FILE_S1='shift.dat',
OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,
cal_lm={int(observables)},cal_cf={int(observables)},NO_two_cf={len(pairs)},FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-14,minitr={minitr},maxitr=1000,itrint=1
/
&input_TRLan
NOE=2,NOK=4,NOM=8,maxitr=1000,lnc_ene_conv=1d-13,i_vec_min=1,i_vec_max=2
/
'''
    (work/'input.dat').write_text(text)
    result = {'name': name, 'algorithm': alg, 'threads': threads, 'cache': cache,
              'expected_dimension': len(values), 'expected_eigenvalues': values.tolist(),
              'input_sha256': {f.name: digest(f) for f in work.glob('*.dat')}}
    try:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        with (work/'input.dat').open('rb') as inp, (work/'stdout.log').open('wb') as log:
            completed = subprocess.run([str(exe)], cwd=work, env=env, stdin=inp,
                                       stdout=log, stderr=subprocess.STDOUT, timeout=120)
        result['exit_code'] = completed.returncode
        assert completed.returncode == 0, f'exit {completed.returncode}'
        out = (work/'stdout.log').read_text()
        assert not re.search(r'\b(?:nan|infinity)\b', out, re.I), 'nonfinite numerical diagnostic in stdout'
        dim = int(re.search(r'THS\(k\)\s*=\s*(\d+)', out).group(1))
        result['dimension'] = dim
        assert dim == len(values), f'sector dimension {dim} != {len(values)}'
        actual = np.loadtxt(work/'eigenvalues.dat', ndmin=2)[:, 1]
        assert np.isfinite(actual).all(), 'nonfinite eigenvalues'
        result['eigenvalues'] = actual.tolist()
        expected = values if alg == 3 else values[:(2 if alg == 2 else 1)]
        assert actual.shape == expected.shape, 'incorrect eigenvalue count'
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=RTOL)
        if observables:
            lm = np.loadtxt(work/'local_mag.dat', ndmin=2)[:, 1:4]
            cf = np.loadtxt(work/'two_body_cf_xyz.dat', ndmin=2)[:, 2:11]
            expected_lm = np.array([[np.vdot(ground, a @ ground).real for a in s] for s in spin])
            expected_cf = np.array([[np.vdot(ground, spin[i][a] @ spin[j][b] @ ground).real
                                     for b in range(3) for a in range(3)] for i,j in pairs])
            # Text output order is xx,yx,zx,xy,yy,zy,xz,yz,zz.
            np.testing.assert_allclose(lm, expected_lm, atol=ATOL, rtol=RTOL)
            np.testing.assert_allclose(cf, expected_cf, atol=ATOL, rtol=RTOL)
            result['local_magnetization'] = lm.tolist()
            result['correlations_xyz'] = cf.tolist()
        result['status'] = 'passed'
    except (AssertionError, OSError, ValueError, AttributeError, subprocess.TimeoutExpired) as error:
        result.update(status='failed', error=str(error))
    print(f"{result['status'].upper()}: {name}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    exe = args.exe.resolve(strict=True)
    c = [1., .8, .7, .2, .1, .3, .1, .3, -.2]
    fields = [[.13*(i+1), -.17+.03*i, .21-.07*i] for i in range(4)]
    bonds = [(0,1,c),(1,2,[x*.7 for x in c]),(2,3,[x*.4 for x in c])]
    cases = []
    for threads in (1,2):
        for alg, cache in [(3,-1),(1,-1),(1,0),(1,3),(2,-1)]:
            cases.append(run_case(exe,args.output,f'xyz-alg{alg}-cache{cache}-omp{threads}',
                                  4,bonds,fields,alg=alg,cache=cache,threads=threads))
    heisenberg = [1.,1.,1.,0.,0.,0.,0.,0.,0.]
    cases.append(run_case(exe,args.output,'early-check-lanczos',4,bonds,fields,alg=1,minitr=1))
    cases.append(run_case(exe,args.output,'singlet-dimer',2,[(0,1,heisenberg)],[[0,0,0]]*2))
    for alg in (3,1,2):
        cases.append(run_case(exe,args.output,f'zero-alg{alg}',4,[],[[0,0,0]]*4,
                              alg=alg,observables=False))
    cases.append(run_case(exe,args.output,'vacuum-lanczos',2,[],[[0,0,1]]*2,
                          alg=1,hi=0,observables=False))
    ring = [(i,(i+1)%6,heisenberg) for i in range(6)]
    for momentum in (1,3):
        for threads in (1,2):
            cases.append(run_case(exe,args.output,f'ring-m{momentum}-omp{threads}',6,ring,
                                  [[0,0,.2]]*6,lo=2,hi=2,momentum=momentum,
                                  threads=threads,observables=False))
    # Nontrivial fixed-excitation sector without a symmetry projection.
    cases.append(run_case(exe,args.output,'truncated-dense',4,bonds,fields,lo=1,hi=2))
    report = {'schema':1, 'numpy':np.__version__, 'executable_sha256':digest(exe),
              'tolerances':{'absolute':ATOL,'relative':RTOL}, 'cases':cases,
              'passed':sum(c['status']=='passed' for c in cases),
              'failed':sum(c['status']=='failed' for c in cases)}
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(f"{report['passed']} passed; {report['failed']} failed; {len(cases)} executed")
    return int(report['failed'] != 0)


if __name__ == '__main__':
    sys.exit(main())
