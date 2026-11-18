#!/usr/bin/env python3
"""Negative and positive source-level guards for noncommuting orbit acceleration.

Real one-dimensional little-group characters are checked against independently
computed Burnside dimensions. Invalid group normal forms/characters must fail
closed to the pre-existing exhaustive basis builder; they must never reach the
new atomic-orbit path. Python standard library only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile

from cpc_openmp_orbit_regression import compile_probe, expected_dimension, run


def square_perm(side, operation):
    def f(i):
        x,y=i%side,i//side
        if operation=='tx':return (x+1)%side+side*y
        if operation=='ty':return x+side*((y+1)%side)
        if operation=='r90':return (-y)%side+side*x
        if operation=='mirror':return x+side*((-y)%side)
        raise ValueError(operation)
    return [f(i) for i in range(side*side)]


def controls():
    n=12
    t=[(i+1)%n for i in range(n)]
    reflect=[(-i)%n for i in range(n)]
    for k, parity in ((0,0),(n//2,1)):
        yield (f'ring12-valid-k{k}-p{parity}', [1]*n, [.5]*n, 2, 2,
               [(t,n,k),(reflect,2,parity)], True)
    yield ('ring12-invalid-inversion-k1', [1]*n, [.5]*n, 2, 2,
           [(t,n,1),(reflect,2,0)], False)
    # X momentum is not C4-invariant; the C4v character is inconsistent.
    side=4
    yield ('square4-invalid-X-C4v', [1]*(side*side), [.5]*(side*side), 3, 3,
           [(square_perm(side,'tx'),side,side//2),
            (square_perm(side,'ty'),side,0),
            (square_perm(side,'r90'),4,0),
            (square_perm(side,'mirror'),2,0)], False)
    # 4-cycle plus a transposition generates S4 (order 24); the proposed
    # cyclic normal form has only 8 maps and is not closed under composition.
    cycle=[1,2,3,0]
    swap01=[1,0,2,3]
    yield ('four-site-nonclosed-S4', [1]*4, [.5]*4, 2, 2,
           [(cycle,4,0),(swap01,2,0)], False)
    yield ('ring12-duplicate-noncommuting-normal-form', [1]*n, [.5]*n, 2, 2,
           [(t,n,0),(reflect,2,0),([(i+2)%n for i in range(n)],n//2,0)],False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--threads',type=int,default=1)
    args=p.parse_args()
    source=args.source.resolve(strict=True)
    root=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='ed2-noncommuting-guard-'))
    if args.output:root.mkdir(parents=True,exist_ok=False)
    probe=compile_probe(source,root/'build',args.compiler,checked=True).resolve(strict=True)
    rows=[]
    for label,caps,spins,lo,hi,gens,valid in controls():
        case=(label,caps,spins,lo,hi,gens)
        actual=run(probe,root/label,case,args.threads)
        fast=actual is not None and actual['parallel_path']
        assert fast==valid, (label,'fast path',fast,'expected',valid)
        if valid:
            expected=expected_dimension(caps,lo,hi,gens)
            assert actual['dim']==expected,(label,actual['dim'],expected)
        else:
            expected=None # Invalid character/group: no Burnside sector claim.
        rows.append(dict(label=label,valid=valid,atomic_fast_path=fast,
                         returned_dimension=actual['dim'] if actual else None,
                         independent_dimension=expected))
        print('PASS',label,'accepted' if fast else 'safe fallback',flush=True)
    assert len(rows)==6
    report=dict(schema=1,passed=True,source_sha256=hashlib.sha256(
                (source/'source'/'state_lists.f90').read_bytes()).hexdigest(),
                checks=rows,caution='A fallback is not proof that invalid symmetry input defines an invariant sector')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('NONCOMMUTING GROUP GUARD 6/6 PASS',flush=True)


if __name__=='__main__':
    main()
