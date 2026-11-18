#!/usr/bin/env python3
"""Independent finite-group and physical model oracle for ED2 simple cubic.

Three lattice translations and three signed-axis point generators give a
unique normal form of all 48 O_h point operations using ED2's six slots:
A = inversion*C4x (order 4), B = C4z^-1 (order 4), C = C3[111] (order 3).
Only Gamma-point one-dimensional O_h characters are encoded.
No scipy/numpy or scheduler needed. Python standard library only.
"""
import argparse
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path

IDENTITY = (1, 2, 3)
IMPROPER_X4 = (-1, -3, 2)   # (-x,z,-y)
PROPER_Z4 = (-2, 1, 3)     # (y,-x,z)
PROPER_C3 = (2, 3, 1)      # (z,x,y)
SECTORS = ("gamma-t", "oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u")
CHARACTERS = {
    "oh-a1g": (0, 0, 0),
    "oh-a2g": (2, 2, 0),
    "oh-a1u": (2, 0, 0),
    "oh-a2u": (0, 2, 0),
}


def apply_signed(mat, vector):
    result = [0, 0, 0]
    for input_axis, axis in enumerate(mat):
        result[abs(axis) - 1] = (1 if axis > 0 else -1) * vector[input_axis]
    return tuple(result)


def mul(mat1, mat2):
    """Matrix composition mat1 after mat2."""
    return tuple((1 if i > 0 else -1) * mat1[abs(i) - 1] for i in mat2)


def powers(mat, order):
    entries = [IDENTITY]
    for _ in range(1, order):
        entries.append(mul(mat, entries[-1]))
    assert mul(mat, entries[-1]) == IDENTITY
    return entries


@lru_cache(maxsize=1)
def full_oh():
    aa, bb, cc = (powers(mat, order) for mat, order in
                  ((IMPROPER_X4, 4), (PROPER_Z4, 4), (PROPER_C3, 3)))
    mapping = {}
    for a, b, c in itertools.product(range(4), range(4), range(3)):
        mat = mul(mul(aa[a], bb[b]), cc[c])
        assert mat not in mapping, "O_h normal form has duplicate operations"
        mapping[mat] = (a, b, c)
    assert len(mapping) == 48
    permutations = {tuple(s * p for s, p in zip(signs, perm))
                    for perm in itertools.permutations((1, 2, 3))
                    for signs in itertools.product((-1, 1), repeat=3)}
    assert set(mapping) == permutations, "O_h normal form is not surjective"
    # Independently enforce that generator phase assignments are real
    # homomorphisms for all pairs, not merely empirical sector dimensions.
    for mode, indices in CHARACTERS.items():
        def character(mat):
            return (-1) ** sum((n // 2) * e for n, e in
                               zip(indices[:2], mapping[mat][:2]))
        assert all(character(mul(a, b)) == character(a) * character(b)
                   for a in mapping for b in mapping), mode
        assert character(IDENTITY) == 1
    return mapping


def index(x, y, z, side):
    return (x % side) + side * ((y % side) + side * (z % side))


def site_xyz(site, side):
    return (site % side, site // side % side, site // (side * side))


def site_permutation(side, transform):
    return tuple(index(*transform(*site_xyz(i, side)), side)
                 for i in range(side ** 3))


def symmetry_generators(side, symmetry):
    if side < 4 or side % 2:
        raise ValueError("Even cubic side >=4 is required for this benchmark")
    if symmetry not in SECTORS:
        raise ValueError(f"Invalid cubic symmetry {symmetry}")
    tx = site_permutation(side, lambda x, y, z: (x + 1, y, z))
    ty = site_permutation(side, lambda x, y, z: (x, y + 1, z))
    tz = site_permutation(side, lambda x, y, z: (x, y, z + 1))
    gens = [(tx, side, 0, "Tx"), (ty, side, 0, "Ty"),
            (tz, side, 0, "Tz")]
    if symmetry != "gamma-t":
        ma, mb, mc = CHARACTERS[symmetry]
        for mat, order, phase, name in (
            (IMPROPER_X4, 4, ma, "I*C4x"),
            (PROPER_Z4, 4, mb, "C4z^-1"),
            (PROPER_C3, 3, mc, "C3[111]"),
        ):
            gens.append((site_permutation(side,
                         lambda x, y, z, m=mat: apply_signed(m, (x, y, z))),
                         order, phase, name))
    return gens


def canonical_bonds(side):
    # One positive axial direction per axis and six face-diagonal directions
    # account for all O_h-symmetric nearest and next-nearest pairs once.
    axial = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    faces = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1),
             (0, 1, 1), (0, 1, -1)]
    bonds = []
    for z in range(side):
        for y in range(side):
            for x in range(side):
                i = index(x, y, z, side)
                for dx, dy, dz in axial:
                    bonds.append((i, index(x + dx, y + dy, z + dz, side), 1., .77))
                for dx, dy, dz in faces:
                    bonds.append((i, index(x + dx, y + dy, z + dz, side), .31, .19))
    assert len(bonds) == 9 * side ** 3
    assert len({(min(i, j), max(i, j)) for i, j, _, _ in bonds}) == len(bonds)
    return bonds


def check_hamiltonian(side):
    bonds = canonical_bonds(side)
    target = sorted((min(i, j), max(i, j), xy, zz) for i, j, xy, zz in bonds)
    for permutation, _, _, name in symmetry_generators(side, "oh-a1g"):
        moved = sorted((min(permutation[i], permutation[j]),
                        max(permutation[i], permutation[j]), xy, zz)
                       for i, j, xy, zz in bonds)
        assert moved == target, f"Cubic J1-J2 XXZ is not invariant under {name}"
    return True


def add_mod(a, b, side):
    return tuple((i + j) % side for i, j in zip(a, b))


def axis_index(xyz, side):
    return index(*xyz, side)


@lru_cache(maxsize=20)
def sector_dimensions(side, down=3):
    if down != 3:
        raise ValueError("Cycle-based physical oracle currently checks three magnons")
    full_oh()
    points = list(itertools.product(range(side), repeat=3))
    size = side ** 3
    assert size >= 64
    sums = {k: 0 for k in SECTORS}
    for mat, (a, b, c) in full_oh().items():
        freq = []
        power = IDENTITY
        for _ in (1, 2, 3):
            power = mul(mat, power)
            h = [0] * size
            for site in points:
                rotated = apply_signed(power, site)
                displacement = tuple((x - y) % side
                                     for x, y in zip(site, rotated))
                h[axis_index(displacement, side)] += 1
            freq.append(h)
        square = mul(mat, mat)
        for shift in points:
            trans2 = add_mod(shift, apply_signed(mat, shift), side)
            trans3 = add_mod(trans2, apply_signed(square, shift), side)
            f1 = freq[0][axis_index(shift, side)]
            two = freq[1][axis_index(trans2, side)] - f1
            three = freq[2][axis_index(trans3, side)] - f1
            assert two >= 0 and two % 2 == 0 and three >= 0 and three % 3 == 0
            fixed_triples = math.comb(f1, 3) + f1 * (two // 2) + (three // 3)
            sums["oh-a1g"] += fixed_triples
            sums["oh-a2g"] += (-1) ** (a + b) * fixed_triples
            sums["oh-a1u"] += (-1) ** a * fixed_triples
            sums["oh-a2u"] += (-1) ** b * fixed_triples
            if mat == IDENTITY:
                sums["gamma-t"] += fixed_triples
    result = {}
    for name, value in sums.items():
        group_size = size if name == "gamma-t" else size * 48
        assert value % group_size == 0, (side, name, value, group_size)
        result[name] = value // group_size
        assert result[name] > 0
    result["raw_dim"] = math.comb(size, down)
    return result


def write_input(work, side, down, gens, bonds, maxiter=2000, wk_dim=-1):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=False)
    n = side ** 3
    (work / "spins.dat").write_text("0.5\n" * n)
    (work / "caps.dat").write_text("1\n" * n)
    (work / "fields.dat").write_text("")
    (work / "bonds.dat").write_text("".join(
        f"{i+1} {j+1} {xy:g} {xy:g} {zz:g} 0 0 0 0 0 0\n"
        for i, j, xy, zz in bonds))
    lines = [f"NOS={n},NODmin={down},NODmax={down},ALG=1,",
             f"NO_one=0,NO_two={len(bonds)},wk_dim={wk_dim},MNTE=-1,",
             "FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',",
             "FILE_spin='spins.dat',FILE_NODmax='caps.dat',",
             "OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0,"]
    for i in range(6):
        if i < len(gens):
            perm, order, phase, label = gens[i]
            fname = f"symmetry{i+1}.dat"
            (work / fname).write_text("".join(f"{j+1}\n" for j in perm))
            lines.append(f"L{i+1}={order},M{i+1}={phase},FILE_S{i+1}='{fname}',")
        else:
            lines.append(f"L{i+1}=1,M{i+1}=0,")
    (work / "input.dat").write_text(
        "&input_parameters\n" + "\n".join(lines) + "\n/\n"
        f"&input_lancz\nlnc_ene_conv=1d-11,minitr=5,maxitr={maxiter},itrint=5\n/\n"
        f"&input_TRLan\nNOE=2,NOK=3,NOM=8,maxitr={maxiter},"
        "lnc_ene_conv=1d-11,i_vec_min=1,i_vec_max=1\n/\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sizes', default='4,6,8,10')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    report = dict(schema=1, geometry='periodic simple-cubic 3D',
                  site_symmetry='O_h, all 48 point operations',
                  representation_scope='four one-dimensional irreps only',
                  results=[])
    for s in (int(x) for x in args.sizes.split(',')):
        assert s >= 4 and s % 2 == 0
        assert check_hamiltonian(s)
        data = sector_dimensions(s)
        print(f'PASS: {s}x{s}x{s} sites={s**3}, dimension={data}', flush=True)
        report['results'].append(dict(side=s, sites=s**3, bond_count=9*s**3,
                                      dimensions=data, checked_bond_invariance=True))
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
