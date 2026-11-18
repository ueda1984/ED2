#!/usr/bin/env python3
"""Verify local-bond sparse cubic ED2 cache against every original cache entry.

Original and local sparse paths run the IDENTICAL anchored cubic basis,
same point-group certificate, representative resolver, input deck and spins.
Require all offdiagonal row indices and complex-element IEEE-754 bits match
*exactly in legacy input-bond order*. Allow only bounded diagonal rounding
from summing the all-up energy first (max abs difference <= 1e-8).
Check exact projected dims, input SHA, energies and certified Ritz residuals.
This is an intentionally stronger check than comparing only ground energies.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import tempfile

from cpc_cubic_pointgroup import (canonical_bonds,sector_dimensions,
                                  symmetry_generators,write_input)
from cpc_pointgroup_profile import run_case


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576),b''):
            h.update(chunk)
    return h.hexdigest()


def frombits(hexdigits):
    return struct.unpack('>d',struct.pack('>Q',int(hexdigits,16)))[0]


def compare_caches(old_path,new_path,dimension,diagonal_max=1e-8):
    compared=0
    offdiag=0
    diag=0
    maximum=0.0
    with old_path.open('rb') as a, new_path.open('rb') as b:
        for left,right in itertools.zip_longest(a,b):
            assert left is not None and right is not None, 'Cache proof line count disagreement'
            x=left.decode('ascii').strip().split()
            y=right.decode('ascii').strip().split()
            assert len(x)==len(y)==5, (compared,x,y)
            assert x[:3]==y[:3], ('Cubic local bond index or order mismatch',compared,x,y)
            compared+=1
            if x[1]=='1':
                diag+=1
                diff=abs(frombits(x[3])-frombits(y[3]))
                imag=abs(frombits(x[4])-frombits(y[4]))
                assert imag==0.0
                assert diff<=diagonal_max, ('Diagonal Hamiltonian mismatch',x,y,diff)
                maximum=max(maximum,diff)
            else:
                offdiag+=1
                assert left==right, (
                    'Offdiagonal complex matrix bit/order mismatch',compared,x,y)
    assert diag==dimension, ('Missing diagonal proof',diag,dimension)
    assert offdiag>0
    return dict(rows=diag,offdiag=offdiag,
                maximum_diagonal_difference=maximum,compared_lines=compared)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--driver',type=Path,required=True)
    p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--output',type=Path)
    p.add_argument('--sides',default='4,6,10')
    p.add_argument('--sectors',default='oh-a1g,oh-a2g,oh-a1u,oh-a2u')
    p.add_argument('--threads',default='1')
    p.add_argument('--timeout',type=int,default=900)
    p.add_argument('--skip-lanczos',action='store_true')
    args=p.parse_args()
    driver=args.driver.resolve(strict=True)
    exe=args.exe.resolve(strict=True)
    root=args.output or Path(tempfile.mkdtemp(prefix='ed2-cubic-local-cache-'))
    if args.output:
        root.mkdir(parents=True,exist_ok=False)
    sectors=args.sectors.split(',')
    assert sectors and set(sectors)<={'oh-a1g','oh-a2g','oh-a1u','oh-a2u'}
    report=[]
    for side in [int(s) for s in args.sides.split(',')]:
        bonds=canonical_bonds(side)
        dimensions=sector_dimensions(side)
        for sector in sectors:
            dim=dimensions[sector]
            for t in (int(x) for x in args.threads.split(',')):
                pair=[]
                for fast in (False,True):
                    work=root/'s{}-{}-t{}-{}'.format(
                        side,sector,t,'local' if fast else 'original')
                    write_input(work,side,3,symmetry_generators(side,sector),
                                bonds,mnte=64)
                    env=dict(os.environ)
                    env.update(ED2_CUBIC_FAST_GROUP_VALIDATION='1',
                               ED2_CUBIC_ANCHORED_BASIS='1',
                               ED2_CUBIC_ANCHORED_REP='1',
                               ED2_CUBIC_CACHE_PROOF_DUMP='1',
                               OMP_NUM_THREADS=str(t),OMP_DYNAMIC='FALSE',
                               OPENBLAS_NUM_THREADS='1')
                    env.pop('ED2_FORCE_LEGACY_ORBITS',None)
                    if fast:
                        env['ED2_CUBIC_LOCAL_BOND_CACHE']='1'
                    else:
                        env.pop('ED2_CUBIC_LOCAL_BOND_CACHE',None)
                    with (work/'input.dat').open('rb') as inp, \
                         (work/'cache_driver_stdout.log').open('wb') as output:
                        cp=subprocess.run([str(driver)],stdin=inp,stdout=output,
                                          stderr=subprocess.STDOUT,cwd=work,
                                          env=env,timeout=args.timeout,check=False)
                    log=(work/'cache_driver_stdout.log').read_text(errors='replace')
                    assert cp.returncode==0,(work,cp.returncode,log[-2000:])
                    assert ('Cubic cache builder: local-bond-sparse' in log) is fast
                    assert 'Representative lookup: cubic-anchored-orbits' in log
                    assert 'Basis builder: anchored-cubic' in log
                    assert 'Cubic cache proof dump: projected rows=' in log
                    assert 'Optimal MNTE =' in log
                    dims=re.findall(r'ED2_CACHE_DRIVER_DIM=(\d+)',log)
                    assert len(dims)==1 and int(dims[0])==dim
                    pair.append((digest(work/'input.dat'),
                                 work/'cubic_cache_proof_dump.txt'))
                assert pair[0][0]==pair[1][0], 'Different source input decks'
                check=compare_caches(pair[0][1],pair[1][1],dim)
                row=dict(side=side,sector=sector,threads=t,dimension=dim,**check)
                report.append(row)
                print('CUBIC_LOCAL_CACHE_BITWISE_PASS side={} sector={} t={} D={} '
                      'offdiag={} max_diag_delta={:.3e}'.format(
                          side,sector,t,dim,check['offdiag'],
                          check['maximum_diagonal_difference']),flush=True)
                if side in (4,6) and t==1 and not args.skip_lanczos:
                    observations=[]
                    for fast in (False,True):
                        work=root/'physical-s{}-{}-{}'.format(
                            side,sector,'local' if fast else 'original')
                        write_input(work,side,3,symmetry_generators(side,sector),
                                    bonds,mnte=64)
                        envnames=('ED2_CUBIC_FAST_GROUP_VALIDATION',
                                  'ED2_CUBIC_ANCHORED_BASIS',
                                  'ED2_CUBIC_ANCHORED_REP',
                                  'ED2_CUBIC_LOCAL_BOND_CACHE')
                        prior={k:os.environ.get(k) for k in envnames}
                        try:
                            os.environ['ED2_CUBIC_FAST_GROUP_VALIDATION']='1'
                            os.environ['ED2_CUBIC_ANCHORED_BASIS']='1'
                            os.environ['ED2_CUBIC_ANCHORED_REP']='1'
                            if fast:
                                os.environ['ED2_CUBIC_LOCAL_BOND_CACHE']='1'
                            else:
                                os.environ.pop('ED2_CUBIC_LOCAL_BOND_CACHE',None)
                            observations.append(run_case(exe,work,side**3,3,sector,
                                                         dim,1,'two-stage',args.timeout))
                        finally:
                            for k,v in prior.items():
                                if v is None:
                                    os.environ.pop(k,None)
                                else:
                                    os.environ[k]=v
                    a,b=observations
                    assert a['input_sha256']==b['input_sha256']
                    assert abs(a['energy']-b['energy'])<1e-8
                    assert a['residual']<=1.01e-8 and b['residual']<=1.01e-8
                    print('CUBIC_LOCAL_CACHE_RITZ_PASS side={} sector={} '
                          'delta_E={:.3e}'.format(
                              side,sector,b['energy']-a['energy']),flush=True)
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('CUBIC_LOCAL_CACHE_REGRESSION_PASS {} exact offdiag and bounded diagonal pairs'.format(
        len(report)),flush=True)


if __name__=='__main__':
    main()
