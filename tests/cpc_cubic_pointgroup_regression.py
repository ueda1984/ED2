#!/usr/bin/env python3
"""Independent orbital Hamiltonian/character projection for 4x4x4 cubic ED2.

The dense reference is assembled in the (small) one-dimensional Oh sector,
not the 41,664-dimensional raw space. Ground-state Ritz pairs are checked
against an explicit spin-bitmask XXZ Hamiltonian, not ED2's own routines.
"""
import argparse
from collections import defaultdict
import itertools
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

from cpc_cubic_pointgroup import (SECTORS, full_oh, canonical_bonds,
                                  sector_dimensions, symmetry_generators, write_input)
from cpc_pointgroup_profile import run_case


def site_index(x, y, z, side):
    return x % side + side * ((y % side) + side * (z % side))


def action_matrix(mat, xyz):
    result = [0, 0, 0]
    for axis, value in enumerate(mat):
        result[abs(value) - 1] = (1 if value > 0 else -1) * xyz[axis]
    return result


def group_actions(side):
    """Independently apply affine signed permutations to each site, Gamma only."""
    for point, (a, b, c) in full_oh().items():
        for tz, ty, tx in itertools.product(range(side), repeat=3):
            action = []
            for site in range(side ** 3):
                xyz = (site % side, site // side % side, site // side ** 2)
                u, v, w = action_matrix(point, xyz)
                action.append(site_index(u + tx, v + ty, w + tz, side))
            yield tuple(action), (a, b, c)


def project_orbit_basis(side):
    assert side == 4
    ops = list(group_actions(side))
    assert len(ops) == side ** 3 * 48
    assert len({g for g, _ in ops}) == len(ops), 'Invalid group actions'
    chars = {'oh-a1g': (0, 0), 'oh-a2g': (1, 1),
             'oh-a1u': (1, 0), 'oh-a2u': (0, 1)}
    states = [sum(1 << i for i in c)
              for c in itertools.combinations(range(side ** 3), 3)]
    visited = set()
    projected = {name: [] for name in chars}
    for rep in states:
        if rep in visited:
            continue
        sites = [i for i in range(side**3) if rep >> i & 1]
        assert len(sites) == 3
        counts = {name: defaultdict(int) for name in chars}
        for perm, (a, b, _) in ops:
            mask = sum(1 << perm[i] for i in sites)
            for name, (ca, cb) in chars.items():
                counts[name][mask] += (-1) ** (ca * a + cb * b)
        visited.update(counts['oh-a1g'])
        for name in chars:
            norm = math.sqrt(sum(x*x for x in counts[name].values()))
            if norm:
                projected[name].append({mask: amp/norm for mask, amp in
                                         counts[name].items() if amp})
    target = sector_dimensions(side)
    for name in chars:
        assert len(projected[name]) == target[name]
    assert len(visited) == len(states)
    return projected


def hamiltonian_reference(side, projected):
    bonds = canonical_bonds(side)
    adjacency = [[] for _ in range(side ** 3)]
    for i, j, xy, zz in bonds:
        adjacency[i].append((j, xy, zz))
        adjacency[j].append((i, xy, zz))
    vacuum = sum(zz*0.25 for _, _, _, zz in bonds)
    energies = {}
    for name, cols in projected.items():
        n = len(cols)
        m = np.zeros((n,n))
        for j, col in enumerate(cols):
            action = defaultdict(float)
            for mask, amplitude in col.items():
                occupied = [s for s in range(side ** 3) if mask >> s & 1]
                assert len(occupied) == 3
                delta = 0.
                for i in occupied:
                    for dest, xy, zz in adjacency[i]:
                        if not (mask >> dest & 1):
                            delta -= zz*.5
                            action[mask ^ (1 << i) ^ (1 << dest)] += xy*.5*amplitude
                action[mask] += (vacuum + delta)*amplitude
            for i, vec in enumerate(cols):
                m[i,j] = sum(action.get(mask, 0.)*weight for mask, weight in vec.items())
            norm2 = sum(x*x for x in action.values())
            projected_norm2 = sum(m[:,j]**2)
            assert abs(norm2 - projected_norm2) < 1e-6, ('Sector leakage', name, j, norm2 - projected_norm2)
        assert np.linalg.norm(m-m.T) < 1e-7
        eigen = np.linalg.eigvalsh(m).tolist()
        energies[name] = eigen
    return energies


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',required=True,type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--threads',type=int,default=1)
    args=p.parse_args()
    root=args.output or Path(tempfile.mkdtemp(prefix='ed2-cubic-oh-reference-'))
    root.mkdir(parents=True,exist_ok=False) if args.output else None
    projected=project_orbit_basis(4)
    spectra=hamiltonian_reference(4,projected)
    bonds=canonical_bonds(4)
    for name in projected:
        dims=sector_dimensions(4)[name]
        gens=symmetry_generators(4,name)
        results={}
        for builder in ('fast','forced-legacy'):
            work=root/(name+'-'+builder)
            write_input(work,4,3,gens,bonds)
            if builder=='forced-legacy':os.environ['ED2_FORCE_LEGACY_ORBITS']='1'
            else:os.environ.pop('ED2_FORCE_LEGACY_ORBITS',None)
            obs=run_case(args.exe.resolve(),work,64,3,name,dims,args.threads,'two-stage',180)
            assert abs(obs['energy']-spectra[name][0]) < 1e-8,(name,builder,obs['energy'],spectra[name][0])
            assert obs['residual'] < 1.1e-8
            results[builder]=obs
        assert abs(results['fast']['energy']-results['forced-legacy']['energy'])<1e-9
        print('PASS',name,'D',dims,'E0',spectra[name][0], 'independent H-sector spectrum',len(spectra[name]),flush=True)
    print('CUBIC OH INDEPENDENT PHYSICAL REGRESSION PASS 4/4',flush=True)

if __name__=='__main__':
    main()
