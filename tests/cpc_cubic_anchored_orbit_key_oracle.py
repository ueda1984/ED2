#!/usr/bin/env python3
"""Independent three-magnon Gamma/O_h *orbit key* proof for a cubic lattice.

This module does not replace the historical ED2 representative. An anchored
key identifies an affine orbit in 3*48 candidate transforms, while ED2's
minimal integer-rank representative can be different. A future fast native
lookup must map each key to the original representative and apply the
ratio of one-dimensional characters; using the key as the basis rank is wrong.

Exhaustive 4^3-spin-orbit tests require only Python 3.6+ stdlib.
"""
import argparse
import itertools
import random

from cpc_cubic_pointgroup import (apply_signed, full_oh, index, site_xyz,
                                  sector_dimensions)


def point_operations(side):
    """Signed-permutation site maps and three nontrivial real O_h chars."""
    assert side >= 4 and side % 2 == 0
    result = []
    for mat, (a, b, c) in full_oh().items():
        mapping = tuple(index(*apply_signed(mat, site_xyz(i, side)), side)
                        for i in range(side ** 3))
        result.append((mapping, ((-1) ** (a + b), (-1) ** a, (-1) ** b)))
    assert len(result) == 48
    return result


def anchor_key_and_characters(state, side, operations, xyz):
    """Smallest zero-anchored triple and all chars reaching that exact key.

    Characters for A2g, A1u and A2u are returned as triples of +/-1.
    If the projected orbit is nonzero in one irrep, the corresponding
    character has a unique value at the minimum key for each source state.
    All operations are real signed coordinate permutations. Gamma makes
    the anchoring translations have character +1.
    """
    best = (side ** 3,) * 3
    signs = set()
    for site_map, char in operations:
        positions = [xyz[site_map[i]] for i in state]
        for anchor in positions:
            candidate = tuple(sorted(index((p[0] - anchor[0]) % side,
                                           (p[1] - anchor[1]) % side,
                                           (p[2] - anchor[2]) % side, side)
                                     for p in positions))
            if candidate < best:
                best = candidate
                signs = {char}
            elif candidate == best:
                signs.add(char)
    assert best[0] == 0 and signs
    return best, signs


def historical_full_orbit_minimum(state, side, operations, xyz):
    """Small exact reference: legacy compares sorted sites from the end."""
    best = (side ** 3,) * 3
    for site_map, _ in operations:
        points = [xyz[site_map[i]] for i in state]
        for tx in range(side):
            for ty in range(side):
                for tz in range(side):
                    ranks = tuple(sorted((index(x + tx, y + ty, z + tz, side)
                                          for x, y, z in points), reverse=True))
                    if ranks < best:
                        best = ranks
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--side", type=int, default=4)
    p.add_argument("--samples", type=int, default=250,
                   help="Equivariance transformations in addition to exhaustive orbit census")
    args = p.parse_args()
    side = args.side
    assert side == 4, "Exhaustive proof currently fixed to 4x4x4"
    assert args.samples > 0
    n = side ** 3
    xyz = tuple(site_xyz(i, side) for i in range(n))
    operations = point_operations(side)

    # Counterexample: orbit minimum under ED2's historical rank convention
    # need not be anchored at the origin even when the orbit key is.
    example = (8, 17, 32)
    anchor, _ = anchor_key_and_characters(example, side, operations, xyz)
    original = historical_full_orbit_minimum(example, side, operations, xyz)
    assert original == (16, 9, 3), original
    assert anchor == (0, 5, 34), anchor
    assert tuple(sorted(original)) != anchor
    print("ANCHORED_NOT_HISTORICAL_REP_PASS", example, original, anchor, flush=True)

    # Every one of the 41,664 full basis triples must have exactly one
    # anchored orbit key. Stabilizer-character conflicts are checked
    # *within each state's mapping to its key*, not across states:
    # different members of an orbit can legitimately carry opposite phases.
    orbits = {}
    for state in itertools.combinations(range(n), 3):
        key, chars = anchor_key_and_characters(state, side, operations, xyz)
        allowed = [len({v[j] for v in chars}) == 1 for j in range(3)]
        if key not in orbits:
            orbits[key] = [0, [True, True, True]]
        orbits[key][0] += 1
        orbits[key][1] = [a and b for a, b in zip(orbits[key][1], allowed)]

    known = sector_dimensions(side)
    measured = [sum(record[1][j] for record in orbits.values()) for j in range(3)]
    assert len(orbits) == known["oh-a1g"]
    assert measured == [known["oh-a2g"], known["oh-a1u"], known["oh-a2u"]]
    assert sum(record[0] for record in orbits.values()) == known["raw_dim"]
    print("ANCHORED_EXHAUSTIVE_ORBITS_PASS raw={raw} A1g={a} A2g={b} A1u={c} A2u={d}".format(
        raw=known["raw_dim"], a=len(orbits), b=measured[0], c=measured[1],
        d=measured[2]), flush=True)

    # For each allowed irrep the signs transform equivariantly under O_h.
    # Character is ±1 so inverse character equals the character itself.
    rng = random.Random(1984)
    for _ in range(args.samples):
        state = tuple(sorted(rng.sample(range(n), 3)))
        key, signs = anchor_key_and_characters(state, side, operations, xyz)
        site_map, character = rng.choice(operations)
        shift = [rng.randrange(side) for _ in range(3)]
        moved = tuple(sorted(index(x + shift[0], y + shift[1], z + shift[2], side)
                             for x, y, z in (xyz[site_map[i]] for i in state)))
        moved_key, moved_signs = anchor_key_and_characters(moved, side, operations, xyz)
        assert moved_key == key
        for irrep in range(3):
            left = {v[irrep] for v in signs}
            right = {v[irrep] for v in moved_signs}
            assert len(left) == len(right)
            if len(left) == 1:
                assert next(iter(left)) == character[irrep] * next(iter(right))
    print("ANCHORED_CHARACTER_EQUIVARIANCE_PASS samples={}".format(args.samples),
          flush=True)

    # Every historical representative and an arbitrary translated/rotated
    # orbit member have the same key, but the stored ED2 basis representative
    # must remain the historical one.
    for _ in range(40):
        state = tuple(sorted(rng.sample(range(n), 3)))
        historical = historical_full_orbit_minimum(state, side, operations, xyz)
        original_key, _ = anchor_key_and_characters(state, side, operations, xyz)
        rep_key, _ = anchor_key_and_characters(tuple(sorted(historical)),
                                                side, operations, xyz)
        assert original_key == rep_key
    print("ANCHORED_HISTORICAL_LOOKUP_KEY_PASS samples=40", flush=True)


if __name__ == "__main__":
    main()
