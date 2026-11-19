#!/usr/bin/env python3
"""Factor exact simple-cubic Oh into 4,3,2,2 symmetry generators.

For a periodic cubic lattice, each rotation is a signed coordinate
permutation (x,y,z)->(s0*r[p0],s1*r[p1],s2*r[p2]).
The ordered product C4*C3*C2*I has exactly 48 unique elements,
and equals the full Oh point group. Combined with L^3 translations it
covers the full cubic space group. This factorization lets an external
basis API use *four* named Oh maps rather than 48 redundant generators.

The proof is independent of QuSpin or ED2 and uses exact integer matrices.
"""
import argparse
import itertools
import json
import math


def mat(rows):
    return tuple(tuple(int(x) for x in row) for row in rows)


IDENTITY = mat(((1,0,0),(0,1,0),(0,0,1)))
C4 = mat(((1,0,0),(0,0,-1),(0,1,0)))
C3 = mat(((0,-1,0),(0,0,-1),(1,0,0)))
C2 = mat(((-1,0,0),(0,0,1),(0,1,0)))
INV = mat(((-1,0,0),(0,-1,0),(0,0,-1)))
FACTORS = (("c4block",C4,4),("c3block",C3,3),
           ("c2block",C2,2),("invblock",INV,2))


def multiply(a,b):
    return mat(tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3))
                           for j in range(3)) for i in range(3)))


def power(a,k):
    result=IDENTITY
    for _ in range(k):
        result=multiply(result,a)
    return result


def full_signed_permutation_group():
    result=set()
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1,1),repeat=3):
            result.add(mat(tuple(tuple(signs[i] if j==perm[i] else 0
                                           for j in range(3)) for i in range(3))))
    return result


def factored_point_group():
    result=set()
    for a in range(4):
        for b in range(3):
            for c in range(2):
                for d in range(2):
                    g=IDENTITY
                    for op,k in ((C4,a),(C3,b),(C2,c),(INV,d)):
                        g=multiply(g,power(op,k))
                    result.add(g)
    return result


def verify_factorization():
    full=full_signed_permutation_group()
    factor=factored_point_group()
    if len(full)!=48 or len(factor)!=48 or full!=factor:
        raise AssertionError("4x3x2x2 does not enumerate exactly Oh")
    for _,g,n in FACTORS:
        if power(g,n)!=IDENTITY or any(power(g,k)==IDENTITY
                                      for k in range(1,n)):
            raise AssertionError("incorrect cyclic generator order")
    for a in factor:
        for b in factor:
            if multiply(a,b) not in factor:
                raise AssertionError("Oh group not closed")
    return dict(status="PASS",group="Oh",order=48,
                exact_factor_orders=[4,3,2,2],
                independent_signed_permutations=len(full),
                generator_names=[x[0] for x in FACTORS])


def site_permutation(side, rotation):
    if side<3:
        raise ValueError("cubic side >=3 required")
    n=side**3
    result=[]
    for i in range(n):
        r=(i%side,(i//side)%side,i//(side*side))
        q=tuple(sum(rotation[a][b]*r[b] for b in range(3)) % side
                for a in range(3))
        result.append(q[0]+side*q[1]+side*side*q[2])
    if len(set(result))!=n:
        raise AssertionError("not a bijective site permutation")
    return result


def cubic_oh_blocks(side):
    return {name:(site_permutation(side,matrix),0)
            for name,matrix,_ in FACTORS}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--side",type=int,default=10)
    a=p.parse_args()
    out=verify_factorization()
    if a.self_test:
        for side in (4,10):
            for _,g,k in FACTORS:
                mapping=site_permutation(side,g)
                current=list(range(side**3))
                for _ in range(k):
                    current=[mapping[i] for i in current]
                if current!=list(range(side**3)):
                    raise AssertionError("bad site-map order")
        out["tested_cubic_sides"]=[4,10]
    else:
        out["site_map_length"]=len(next(iter(cubic_oh_blocks(a.side).values()))[0])
    print(json.dumps(out,sort_keys=True))
