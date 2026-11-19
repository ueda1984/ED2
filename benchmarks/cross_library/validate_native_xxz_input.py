#!/usr/bin/env python3
"""Check an ED2 native cubic bond deck against the exact PR16 XXZ model.

This checks ONLY the complete two-spin interaction terms. It is not a proof
of input.dat spin lengths, magnetization window, fields or symmetry sector.
"""
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from export_cubic_xxz import cubic_xxz_bonds, metadata


def read_native(path, n):
    rows=[]
    for lineno,line in enumerate(Path(path).read_text().splitlines(), 1):
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        words=line.split()
        if len(words)!=11:
            raise ValueError("ED2 line %d: expected i j Jx Jy Jz Dx Dy Dz Gx Gy Gz" % lineno)
        try:
            i,j=int(words[0])-1,int(words[1])-1
            val=[float(s) for s in words[2:]]
        except ValueError as e:
            raise ValueError("invalid field on ED2 bond line %d" % lineno) from e
        if not(0 <= i < n and 0 <= j < n and i!=j):
            raise ValueError("ED2 line %d: invalid site IDs" % lineno)
        if not all(math.isfinite(x) for x in val):
            raise ValueError("ED2 line %d: nonfinite coefficient" % lineno)
        if val[0]!=val[1] or any(x!=0.0 for x in val[3:]):
            raise ValueError("ED2 line %d: not real XXZ without DM/Gamma" % lineno)
        rows.append((min(i,j),max(i,j),val[0],val[2]))
    return rows


def validate(side,path):
    expect=Counter((min(i,j),max(i,j),xy,zz)
                   for i,j,xy,zz in cubic_xxz_bonds(side))
    found=Counter(read_native(path,side**3))
    if expect!=found:
        missing=sum((expect-found).values())
        extra=sum((found-expect).values())
        raise ValueError("native bonds differ from PR16: missing=%d extra=%d" %
                         (missing,extra))
    return dict(status="PASS_BONDS_ONLY",bond_count=sum(found.values()),
                model=metadata(side,3)["model"],
                physical_sha256_m3=metadata(side,3)["physical_sha256"],
                note="No assertion yet about spin/field/NOD/symmetry/solver equivalence")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--side",type=int,required=True)
    p.add_argument("--ed2-bonds",required=True)
    a=p.parse_args()
    print(json.dumps(validate(a.side,a.ed2_bonds),sort_keys=True))


if __name__=="__main__":
    main()
