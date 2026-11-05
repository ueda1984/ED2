#!/usr/bin/env python3
"""Independent field-record, solver and full operator validation (test-only NumPy)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import numpy as np
from numerical_baseline import model

ATOL = RTOL = 1e-9
OPERATOR_ATOL = 1e-12
N = 4
COUPLINGS = [1., .8, .7, .2, .1, .3, .1, .3, -.2]
BONDS = [(0,1,COUPLINGS),(1,2,[x*.7 for x in COUPLINGS]),(2,3,[x*.4 for x in COUPLINGS])]
FIELDS = [[.13*(i+1), -.17+.03*i, .21-.07*i] for i in range(N)]
# Sites are zero-based only in this Python reference; ED2 input remains one-based.
FIXTURES = {
    'single-site4': [(3,[.37,-.23,.41])],
    'sparse-ordered': [(0,[.11,-.17,.29]),(2,[-.31,.19,.43])],
    'sparse-permuted': [(2,[-.31,.19,.43]),(0,[.11,-.17,.29])],
    'full-ordered': list(enumerate(FIELDS)),
    'full-permuted': [(i,FIELDS[i]) for i in [2,0,3,1]],
}
CONFIGURATIONS = [(1,c) for c in [0,3,-1]] + [(2,c) for c in [0,3,-1]] + [(3,-1)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reference(records):
    fields = np.zeros((N,3))
    for site, field in records:
        fields[site] += field
    h, spin = model(N,BONDS,fields)
    np.testing.assert_allclose(h,h.conj().T,atol=OPERATOR_ATOL,rtol=0)
    values, vectors = np.linalg.eigh(h)
    require(values[1]-values[0]>1e-6,'fixture ground state must be nondegenerate')
    psi = vectors[:,0]
    lm = np.array([[np.vdot(psi,a@psi).real for a in site] for site in spin])
    cf = np.array([[np.vdot(psi,spin[i][a]@spin[j][b]@psi).real
                    for b in range(3) for a in range(3)]
                   for i in range(N) for j in range(i+1,N)])
    return h, values, lm, cf


def prepare(work, records, alg, cache):
    work.mkdir()
    (work/'bonds.dat').write_text(''.join(f'{i+1} {j+1} '+' '.join(map(str,c))+'\n' for i,j,c in BONDS))
    (work/'fields.dat').write_text(''.join(f'{i+1} '+' '.join(map(str,h))+'\n' for i,h in records))
    (work/'spins.dat').write_text('0.5\n'*N)
    (work/'caps.dat').write_text('1\n'*N)
    (work/'pairs.dat').write_text(''.join(f'{i+1} {j+1}\n' for i in range(N) for j in range(i+1,N)))
    (work/'input.dat').write_text(f'''&input_parameters
NOS={N},NODmin=0,NODmax={N},ALG={alg},
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={len(records)},NO_two={len(BONDS)},wk_dim={cache},MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=1,cal_cf=1,NO_two_cf=6,FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-14,minitr=5,maxitr=1000,itrint=1
/
&input_TRLan
NOE=2,NOK=4,NOM=8,maxitr=1000,lnc_ene_conv=1d-13,i_vec_min=1,i_vec_max=2
/
''')
    return {p.name:sha(p) for p in sorted(work.glob('*.dat'))}


def execute(exe, work, threads):
    env = dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (work/'input.dat').open('rb') as inp, (work/'stdout.log').open('wb') as log:
        completed = subprocess.run([str(exe)],cwd=work,env=env,stdin=inp,
                                   stdout=log,stderr=subprocess.STDOUT,timeout=120)
    return completed.returncode


def compare(actual, expected, label, atol=ATOL, rtol=RTOL):
    require(actual.shape==expected.shape,label+': shape differs')
    require(np.isfinite(actual).all(),label+': nonfinite result')
    np.testing.assert_allclose(actual,expected,atol=atol,rtol=rtol,err_msg=label)
    return float(np.max(np.abs(actual-expected)))


def cli_case(exe, root, name, records, alg, cache, threads):
    work = root/f'{name}-alg{alg}-cache{cache}-omp{threads}'
    result = {'fixture':name,'alg':alg,'cache':cache,'threads':threads,
              'input_sha256':prepare(work,records,alg,cache)}
    _, values, lm, cf = reference(records)
    try:
        result['exit_code'] = execute(exe,work,threads)
        require(result['exit_code']==0,f"exit {result['exit_code']}")
        log = (work/'stdout.log').read_text()
        require(not re.search(r'\b(?:nan|infinity)\b',log,re.I),'nonfinite diagnostic')
        expected = values if alg==3 else values[:(2 if alg==2 else 1)]
        eigenvalues = np.loadtxt(work/'eigenvalues.dat',ndmin=2)
        np.testing.assert_array_equal(eigenvalues[:,0],np.arange(1,len(expected)+1))
        local = np.loadtxt(work/'local_mag.dat',ndmin=2)
        correlations = np.loadtxt(work/'two_body_cf_xyz.dat',ndmin=2)
        np.testing.assert_array_equal(local[:,0],np.arange(1,N+1))
        np.testing.assert_array_equal(correlations[:,:2],[(i+1,j+1) for i in range(N) for j in range(i+1,N)])
        result['errors'] = {'energy':compare(eigenvalues[:,1],expected,'spectrum'),
                            'magnetization':compare(local[:,1:],lm,'magnetization'),
                            'correlation':compare(correlations[:,2:],cf,'correlations')}
        result.update(status='passed',eigenvalues=eigenvalues[:,1].tolist(),
                      magnetization=local[:,1:].tolist(),correlations=correlations[:,2:].tolist())
    except (OSError,ValueError,AssertionError,subprocess.TimeoutExpired) as error:
        result.update(status='failed',error=str(error))
    return result


def matrix_file(path, dim):
    data = np.loadtxt(path,ndmin=2)
    require(data.shape==(dim*dim,4),'operator export shape differs')
    rows, cols = data[:,0].astype(int)-1, data[:,1].astype(int)-1
    require(np.array_equal(data[:,0],rows+1) and np.array_equal(data[:,1],cols+1),'noninteger operator index')
    require(np.all((rows>=0)&(rows<dim)&(cols>=0)&(cols<dim)),'operator index out of range')
    require(len(set(zip(rows,cols)))==dim*dim,'operator indices repeated')
    matrix = np.zeros((dim,dim),complex)
    matrix[rows,cols] = data[:,2]+1j*data[:,3]
    require(np.isfinite(matrix).all(),'nonfinite operator')
    return matrix


def operator_case(exe, root, name, records, threads):
    work = root/f'operator-{name}-omp{threads}'
    result = {'fixture':name,'threads':threads,'input_sha256':prepare(work,records,3,-1)}
    h, _, _, _ = reference(records)
    try:
        result['exit_code'] = execute(exe,work,threads)
        require(result['exit_code']==0,f"exit {result['exit_code']}")
        basis = np.loadtxt(work/'basis.dat',dtype=int,ndmin=2)
        require(basis.shape==(2**N,N+1),'basis shape differs')
        np.testing.assert_array_equal(basis[:,0],np.arange(1,2**N+1))
        permutation = []
        for row in basis[:,1:]:
            occupied = [int(site) for site in row if site!=0]
            require(len(occupied)==len(set(occupied)) and all(1<=i<=N for i in occupied),'invalid occupation')
            permutation.append(sum(1<<(N-site) for site in occupied))
        require(sorted(permutation)==list(range(2**N)),'basis is not complete')
        expected = h[np.ix_(permutation,permutation)]
        dense = matrix_file(work/'dense.dat',2**N)
        errors = {'dense_reference':compare(dense,expected,'dense reference',OPERATOR_ATOL,0),
                  'dense_hermiticity':compare(dense,dense.conj().T,'dense Hermiticity',OPERATOR_ATOL,0)}
        for cache in [0,3,-1]:
            applied = matrix_file(work/f'apply-{cache}.dat',2**N)
            errors[f'cache{cache}_reference'] = compare(applied,expected,'H applied to complete basis',OPERATOR_ATOL,0)
            errors[f'cache{cache}_dense'] = compare(applied,dense,'cached/on-demand versus dense',OPERATOR_ATOL,0)
            errors[f'cache{cache}_hermiticity'] = compare(applied,applied.conj().T,'application Hermiticity',OPERATOR_ATOL,0)
        result.update(status='passed',errors=errors,basis_permutation=permutation)
    except (OSError,ValueError,AssertionError,subprocess.TimeoutExpired) as error:
        result.update(status='failed',error=str(error))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--operator-exe',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    exe = args.exe.resolve(strict=True)
    operator_exe = args.operator_exe.resolve(strict=True)
    if args.output is None:
        out = Path(tempfile.mkdtemp(prefix='ed2-sparse-field-'))
    else:
        out = args.output.resolve()
        out.mkdir(parents=True,exist_ok=False)
    cli, operators, invariances = [], [], []
    for name, records in FIXTURES.items():
        for threads in [1,2]:
            operators.append(operator_case(operator_exe,out,name,records,threads))
            for alg,cache in CONFIGURATIONS:
                result = cli_case(exe,out,name,records,alg,cache,threads)
                cli.append(result)
        print(name,':',sum(c['status']=='passed' for c in cli if c['fixture']==name),'CLI cases passed',flush=True)
    for ordered,permuted in [('sparse-ordered','sparse-permuted'),('full-ordered','full-permuted')]:
        for threads in [1,2]:
            for alg,cache in CONFIGURATIONS:
                pair = [next(c for c in cli if c['fixture']==name and c['threads']==threads and c['alg']==alg and c['cache']==cache) for name in [ordered,permuted]]
                result = {'fixtures':[ordered,permuted],'alg':alg,'cache':cache,'threads':threads}
                try:
                    require(all(c['status']=='passed' for c in pair),'prerequisite solver case failed')
                    result['errors'] = {key:compare(np.array(pair[0][key]),np.array(pair[1][key]),'record permutation '+key)
                                        for key in ['eigenvalues','magnetization','correlations']}
                    result['status'] = 'passed'
                except (ValueError,AssertionError) as error:
                    result.update(status='failed',error=str(error))
                invariances.append(result)
    def counts(cases):
        return {'passed':sum(c['status']=='passed' for c in cases),'failed':sum(c['status']=='failed' for c in cases)}
    report = {'schema':1,'numpy':np.__version__,'executable_sha256':sha(exe),'operator_executable_sha256':sha(operator_exe),
              'harness_sha256':sha(Path(__file__)),'reference_constructor_sha256':sha(Path(__file__).with_name('numerical_baseline.py')),
              'tolerances':{'absolute':ATOL,'relative':RTOL,'operator_absolute':OPERATOR_ATOL},
              'fixtures':{name:{'records':records,'ground_gap':float(np.diff(reference(records)[1][:2])[0])} for name,records in FIXTURES.items()},
              'cli_counts':counts(cli),'operator_counts':counts(operators),'permutation_counts':counts(invariances),
              'cli':cli,'operators':operators,'permutations':invariances}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('CLI:',report['cli_counts'],'operator:',report['operator_counts'],'permutations:',report['permutation_counts'])
    print('Evidence:',out)
    return int(any(report[key]['failed'] for key in ['cli_counts','operator_counts','permutation_counts']))


if __name__=='__main__':
    sys.exit(main())
