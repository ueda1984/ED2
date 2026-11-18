#!/usr/bin/env python3
"""Independent exact colex oracle for axis-pruned cubic three-site translations.

Exhaustively enumerate every distinct three-site configuration at L=4
(C(64,3)=41664), then sample seeded large-L configurations. Both algorithms
are implemented here independently of ED2 Fortran and bitwise proof drivers.
Python >= 3.6; standard library only.
"""
import argparse
import itertools
import random


def delta_axes(sites, side):
    # d[axis][anchor][site]; site and anchor in 0..2.
    return [[[((sites[j][axis] - sites[i][axis]) % side)
              for j in range(3)] for i in range(3)] for axis in range(3)]


def triple(d, side, pi, pj, pk):
    return tuple(sorted(1 + d[0][pi][j] +
                        side * (d[1][pj][j] + side * d[2][pk][j])
                        for j in range(3)))


def canonical_full(d, side):
    best = None
    for pi, pj, pk in itertools.product(range(3), repeat=3):
        cur = triple(d, side, pi, pj, pk)
        if best is None or cur[::-1] < best[::-1]:
            best = cur
    return best


def canonical_pruned(d, side):
    z_high = [max(d[2][pk]) for pk in range(3)]
    z_min = min(z_high)
    best = None
    candidates = 0
    for pk in range(3):
        if z_high[pk] != z_min:
            continue
        y_high = [max(d[1][pj][j] for j in range(3)
                      if d[2][pk][j] == z_min) for pj in range(3)]
        y_min = min(y_high)
        for pj in range(3):
            if y_high[pj] != y_min:
                continue
            x_high = [max(d[0][pi][j] for j in range(3)
                          if d[2][pk][j] == z_min and
                          d[1][pj][j] == y_min) for pi in range(3)]
            x_min = min(x_high)
            for pi in range(3):
                if x_high[pi] != x_min:
                    continue
                cur = triple(d, side, pi, pj, pk)
                if best is None or cur[::-1] < best[::-1]:
                    best = cur
                candidates += 1
    assert 1 <= candidates <= 27
    return best, candidates


def point(site, side):
    return (site % side, (site // side) % side, site // (side * side))


def check(indices, side):
    d = delta_axes([point(index, side) for index in indices], side)
    old = canonical_full(d, side)
    new, count = canonical_pruned(d, side)
    assert old == new, ("Axis-pruned canonicalization changed colex", side,
                        indices, old, new)
    return count


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--random-per-side", type=int, default=15000)
    args = p.parse_args()
    assert args.random_per_side > 0
    counts = []
    for indices in itertools.combinations(range(4 ** 3), 3):
        counts.append(check(indices, 4))
    assert len(counts) == 41664
    print("CUBIC_PRUNED_CANON_EXHAUSTIVE_PASS L=4 cases={} "
          "mean_candidates={:.4f}".format(
              len(counts), sum(counts) / float(len(counts))), flush=True)

    rng = random.Random(1984)
    for side in (6, 8, 10):
        counts = []
        for _ in range(args.random_per_side):
            indices = rng.sample(range(side ** 3), 3)
            counts.append(check(indices, side))
        print("CUBIC_PRUNED_CANON_RANDOM_PASS L={} cases={} "
              "mean_candidates={:.4f} max_candidates={}".format(
                  side, len(counts), sum(counts) / float(len(counts)),
                  max(counts)), flush=True)
    print("CUBIC_PRUNED_CANON_INTEGER_ORACLE_PASS", flush=True)


if __name__ == "__main__":
    main()
