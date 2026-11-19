#!/usr/bin/env python3
"""Small cubic J1-J2 XXZ Γ/Oh A1g eigenpair *independent physical oracle*.

The full spin configuration basis is explicitly enumerated for side=4 only.
A BFS of translation, coordinate reflection and axis-swap generators assigns
every raw 3-magnon configuration to a complete space-group orbit. Each orbit
is normalized as (sum_{config in orbit}|config>)/sqrt(orbit size).
H is independently assembled from local bonds (rather than ED2).
Do not use as a 1000-site performance benchmark.
"""
import argparse
import itertools
import json
import math
import numpy as np
from export_cubic_xxz import cubic_xxz_bonds, metadata
from verify_orbits import generators


def projected_reference(side=4, ndown=3):
    n = side ** 3
    if side != 4 or ndown not in (2, 3) or math.comb(n, ndown) > 50000:
        raise ValueError("small independent reference limited to side=4, m=2/3")
    bonds = cubic_xxz_bonds(side)
    adj = [[] for _ in range(n)]
    jzz_pair = {}
    for a,b,xy,zz in bonds:
        adj[a].append((b,xy,zz))
        adj[b].append((a,xy,zz))
        jzz_pair[tuple(sorted((a,b)))] = zz
    degree_zz = [sum(zz for _,_,zz in nbrs) for nbrs in adj]
    eall = sum(zz for _,_,_,zz in bonds) / 4.0
    where = {}
    sizes = []
    gs = generators((side,side,side))
    for root in itertools.combinations(range(n),ndown):
        if root in where:
            continue
        idx = len(sizes)
        stack = [root]
        where[root] = idx
        length = 0
        while stack:
            current = stack.pop()
            length += 1
            for g in gs:
                nxt = tuple(sorted(g[x] for x in current))
                if nxt not in where:
                    where[nxt] = idx
                    stack.append(nxt)
        sizes.append(length)
    assert len(where) == math.comb(n,ndown)
    assert sum(sizes) == len(where)
    d = len(sizes)
    h = np.zeros((d,d),dtype=np.float64)
    # H|configuration>: diagonal = all-up background - one-down endpoint
    # losses + double-down correction; swap a down to an adjacent up with
    # amplitude J_xy/2. Project with exact orbit normalizations.
    for state, col in where.items():
        state_set = set(state)
        diagonal = eall - 0.5 * sum(degree_zz[a] for a in state)
        diagonal += sum(jzz_pair.get((a,b), 0.0)
                        for a,b in itertools.combinations(state,2))
        h[col,col] += diagonal / sizes[col]
        for a in state:
            for b,xy,_ in adj[a]:
                if b in state_set:
                    continue
                moved = tuple(sorted((state_set-{a}) | {b}))
                row = where[moved]
                h[row,col] += xy * 0.5 / math.sqrt(sizes[col]*sizes[row])
    asymmetry = float(np.max(np.abs(h-h.T)))
    if asymmetry > 1e-10:
        raise AssertionError("projected Hamiltonian not symmetric")
    e, v = np.linalg.eigh(h)
    residual = float(np.linalg.norm(h@v[:,0] - e[0]*v[:,0]))
    return dict(metadata(side, ndown), projected_dimension=d,
                symmetry="Gamma/O_h/A1g from full real-space generator orbit",
                raw_states_enumerated=len(where),
                E0=float(e[0]), E1=float(e[1]),
                matrix_asymmetry_max=asymmetry,
                direct_residual=residual,
                status="PASS" if residual < 1e-10 else "BAD_RESIDUAL")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--ndown",type=int,choices=[2,3],default=3)
    p.add_argument("--self-test",action="store_true")
    a = p.parse_args()
    o = projected_reference(4, a.ndown)
    if a.self_test:
        if a.ndown == 3:
            assert o["projected_dimension"] == 44
            assert abs(o["E0"] - 40.345908939159166) < 1e-9
        assert o["status"] == "PASS"
    print(json.dumps(o,sort_keys=True))
