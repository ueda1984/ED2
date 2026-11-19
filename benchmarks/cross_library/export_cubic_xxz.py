#!/usr/bin/env python3
"""Export the *exact physical* J1-J2 XXZ spin-half cubic model from ED2 PR16.

Source reference: tests/cpc_cubic_pointgroup.py::canonical_bonds,
ueda1984/ED2, Draft PR16 commit 9692d8b1fd903bcb7b1b67a5f180c57430702510.
This exports source-derived parameters but does not claim byte-identical
ED2 input decks or an independently verified inter-library energy.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

AXIAL = ((1,0,0),(0,1,0),(0,0,1))
FACES = ((1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1))
SOURCE = "ueda1984/ED2 PR16 9692d8b1fd903bcb7b1b67a5f180c57430702510"


def cubic_xxz_bonds(side):
    """Preserve native ED2 canonical_bonds positive-axis/face ordering."""
    if side < 4 or side % 2:
        raise ValueError("cubic side must be even, >=4 (PR16 certified scope)")
    def index(x, y, z):
        return (x % side) + side * ((y % side) + side * (z % side))
    bonds = []
    for z in range(side):
        for y in range(side):
            for x in range(side):
                i = index(x, y, z)
                for dx, dy, dz in AXIAL:
                    bonds.append((i, index(x+dx, y+dy, z+dz), 1.0, 0.77))
                for dx, dy, dz in FACES:
                    bonds.append((i, index(x+dx, y+dy, z+dz), 0.31, 0.19))
    assert len(bonds) == 9*side**3
    assert len({(min(i,j), max(i,j)) for i,j,_,_ in bonds}) == len(bonds)
    return bonds


def metadata(side, ndown):
    if ndown not in (2,3):
        raise ValueError("benchmark only supports m=2 or m=3")
    n = side**3
    ordered = cubic_xxz_bonds(side)
    canonical = sorted((min(i,j),max(i,j),jxy,jzz) for i,j,jxy,jzz in ordered)
    data = {"model":"periodic_spin_half_cubic_J1J2_XXZ", "shape":[side]*3,
            "ndown":ndown, "site_order":"x-fastest", "field":0.0,
            "bonds":[list(e) for e in canonical],
            "spin_operator_convention":"SdotS uses spin (not Pauli) normalization"}
    sha = hashlib.sha256(json.dumps(data,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    return {"model":data["model"],"shape":[side]*3,"sites":n,"ndown":ndown,
            "u1_dimension":math.comb(n,ndown),
            "bond_count":len(ordered),"j1_bonds":3*n,"j2_bonds":6*n,
            "j1_xy":1.0,"j1_zz":0.77,"j2_xy":0.31,"j2_zz":0.19,
            "field":0.0,"site_order":"x-fastest",
            "physical_sha256":sha,"source":SOURCE,
            "all_up_energy":(3*n*0.77+6*n*0.19)/4.0}


def export(side, ndown, prefix):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    meta=metadata(side,ndown)
    bonds=cubic_xxz_bonds(side)
    data="".join("%d %d %.17g %.17g\n"%(i,j,jxy,jzz)
                 for i,j,jxy,jzz in bonds)
    bondfile=prefix.with_suffix(".xxz_bonds")
    manifestfile=prefix.with_suffix(".json")
    bondfile.write_text(data,encoding="ascii")
    manifest=dict(meta)
    manifest.update(
        bond_file=bondfile.name,
        bond_file_sha256=hashlib.sha256(data.encode("ascii")).hexdigest(),
        input_columns="0-based endpoint_i endpoint_j Jxy Jzz",
        quspin_terms="zz=Jzz, +-=Jxy/2, -+=Jxy/2; pauli=0",
        qualified_cubic_specialized_pr="Ndown=3, Gamma+Oh 1D irrep, exact 9N bonds",
        safety="export only; no Hilbert space allocated or solver run")
    manifestfile.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return manifest


def self_test():
    for side in (4,6,10):
        b=cubic_xxz_bonds(side)
        m=metadata(side,3)
        assert len(b)==9*side**3
        assert m["j1_bonds"]+m["j2_bonds"]==len(b)
        assert abs(m["all_up_energy"]-.8625*side**3)<1e-9
    assert abs(metadata(10,3)["all_up_energy"]-862.5)<1e-10
    assert metadata(10,3)["u1_dimension"]==166167000
    # No duplicate bonds and correct O_h/translation generator action.
    from verify_orbits import generators
    for side in (4,10):
        bonds=cubic_xxz_bonds(side)
        original={(min(i,j),max(i,j),jxy,jzz) for i,j,jxy,jzz in bonds}
        for g in generators((side,side,side)):
            moved={(min(g[i],g[j]),max(g[i],g[j]),jxy,jzz)
                   for i,j,jxy,jzz in bonds}
            assert original==moved
    print("PASS: ED2 PR16 J1-J2 XXZ canonical cubic bonds, O_h/translation generators, Eallup")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--side",type=int,default=10)
    p.add_argument("--ndown",type=int,default=3)
    p.add_argument("--output-prefix")
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test()
    elif a.output_prefix:
        print(json.dumps(export(a.side,a.ndown,a.output_prefix),sort_keys=True))
    else:
        print(json.dumps(metadata(a.side,a.ndown),sort_keys=True))
