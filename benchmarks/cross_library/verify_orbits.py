#!/usr/bin/env python3
"""Exhaustively enumerate *small* lattice-space-group orbits, independently
of Burnside fixed-point/cycle arithmetic. Intended solely for scientific QA.
"""
import argparse
import itertools
import json
import math
from reference import parse_shape
from symmetry_dimensions import count_blocks


def generators(shape):
    dim = len(shape)
    n = math.prod(shape)
    strides = [math.prod(shape[:a]) for a in range(dim)]
    coords = [tuple((site // strides[a]) % shape[a] for a in range(dim))
              for site in range(n)]
    def encode(c):
        return sum(c[a] * strides[a] for a in range(dim))
    maps = []
    # Translations, individual coordinate sign changes and axis swaps
    # generate Z_L^d semidirect the signed permutation point group.
    for axis in range(dim):
        maps.append(tuple(encode(tuple((v + (1 if a == axis else 0)) % shape[a]
                                       for a, v in enumerate(c))) for c in coords))
        maps.append(tuple(encode(tuple(((-v if a == axis else v) % shape[a])
                                       for a, v in enumerate(c))) for c in coords))
    for a in range(dim):
        for b in range(a+1, dim):
            if shape[a] != shape[b]:
                continue
            tab = []
            for c in coords:
                x = list(c)
                x[a], x[b] = x[b], x[a]
                tab.append(encode(x))
            maps.append(tuple(tab))
    for m in maps:
        assert sorted(m) == list(range(n))
    return maps


def brute_orbits(shape, m, max_states=50000):
    n = math.prod(shape)
    if m not in (2, 3) or math.comb(n, m) > max_states:
        raise ValueError("exhaustive orbit enumeration exceeds cap")
    gs = generators(shape)
    seen = set()
    orbits = 0
    for root in itertools.combinations(range(n), m):
        if root in seen:
            continue
        orbits += 1
        stack = [root]
        seen.add(root)
        while stack:
            current = stack.pop()
            for mapping in gs:
                nxt = tuple(sorted(mapping[s] for s in current))
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
    assert len(seen) == math.comb(n, m)
    return orbits


def verify(shape, m):
    found = brute_orbits(shape, m)
    predicted = count_blocks(shape)["m%d" % m]["gamma_point_A1"]
    if found != predicted:
        raise AssertionError("Burnside %d != brute orbit %d" %
                             (predicted, found))
    return {"shape": list(shape), "ndown": m, "orbits": found,
            "status": "PASS", "method": "explicit generator BFS vs Burnside"}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--shape", default="4,4,4")
    p.add_argument("--ndown", type=int, default=3)
    a = p.parse_args()
    if a.self_test:
        for shape in ((6,), (4, 4), (4, 4, 4)):
            for m in (2, 3):
                print(json.dumps(verify(shape, m), sort_keys=True), flush=True)
        print("PASS: six independently enumerated orbit counts")
    else:
        print(json.dumps(verify(parse_shape(a.shape), a.ndown), sort_keys=True))
