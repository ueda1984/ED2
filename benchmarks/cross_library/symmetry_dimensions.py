#!/usr/bin/env python3
"""Independent Burnside count of spin-half 2/3-magnon space-group A1 blocks.

Point transformations are the signed coordinate permutations that preserve
the lattice shape: C2 in 1D, C4v in 2D square, Oh in 3D cube.
Only the TRIVIAL 1D representation (Gamma / A1) is calculated; it does not
certify other Oh irreps or a Hamiltonian's actual symmetry.
"""
import argparse
import itertools
import json
import math
import time

import numpy as np
from reference import case_metadata, parse_shape


PRIMARY_SHAPES = ((1000,), (32, 32), (10, 10, 10))


def count_blocks(shape, max_sites=2048):
    n = math.prod(shape)
    if n > max_sites:
        raise ValueError("refuse O(N^2) memory above max_sites")
    if len(shape) > 3 or any(s < 3 for s in shape):
        raise ValueError("expect 1D/2D/3D periodic axes >=3")

    t0 = time.perf_counter()
    dims = np.asarray(shape, dtype=np.int32)
    stride = np.asarray([math.prod(shape[:i]) for i in range(len(shape))],
                        dtype=np.int32)
    labels = np.arange(n, dtype=np.int32)
    coords = (labels[:, None] // stride[None, :]) % dims[None, :]
    # shift[t, s] = t+s, using x-fastest labels. O(N^2) only.
    shift = ((coords[None, :, :] + coords[:, None, :]) %
             dims[None, None, :]) @ stride
    gamma_counts = [0, 0]
    a1_counts = [0, 0]
    point_order = 0
    for perm in itertools.permutations(range(len(shape))):
        if any(shape[perm[i]] != shape[i] for i in range(len(shape))):
            continue
        for signs in itertools.product((-1, 1), repeat=len(shape)):
            rotated = ((coords[:, perm] *
                        np.asarray(signs, dtype=np.int32)[None, :]) %
                       dims[None, :]) @ stride
            affine = shift[:, rotated]
            f1 = np.count_nonzero(affine == labels[None, :], axis=1).astype(np.int64)
            g2 = np.take_along_axis(affine, affine, axis=1)
            f2 = np.count_nonzero(g2 == labels[None, :], axis=1).astype(np.int64)
            g3 = np.take_along_axis(affine, g2, axis=1)
            f3 = np.count_nonzero(g3 == labels[None, :], axis=1).astype(np.int64)
            if not (np.all((f2 - f1) % 2 == 0) and
                    np.all((f3 - f1) % 3 == 0)):
                raise AssertionError("invalid cycle decomposition")
            cycle2 = (f2 - f1) // 2
            cycle3 = (f3 - f1) // 3
            # Invariant 2/3-subsets are unions of cycles with sizes 1,2,3.
            fixed2 = f1 * (f1 - 1) // 2 + cycle2
            fixed3 = (f1 * (f1 - 1) * (f1 - 2) // 6 +
                      f1 * cycle2 + cycle3)
            for i, fixed in enumerate((fixed2, fixed3)):
                a1_counts[i] += int(np.sum(fixed))
                if perm == tuple(range(len(shape))) and signs == (1,) * len(shape):
                    gamma_counts[i] += int(np.sum(fixed))
            point_order += 1

    if any(v % n for v in gamma_counts) or any(v % (n * point_order)
                                             for v in a1_counts):
        raise AssertionError("Burnside sum not divisible by group order")
    result = {
        "shape": list(shape),
        "sites": n,
        "point_group_order": point_order,
        "space_group_order": n * point_order,
        "symmetry": "Gamma and trivial (A1) point representation",
        "elapsed_seconds_oracle": round(time.perf_counter() - t0, 5),
    }
    for m, idx in ((2, 0), (3, 1)):
        result["m%d" % m] = {
            "fixed_magnetization": math.comb(n, m),
            "gamma_translation": gamma_counts[idx] // n,
            "gamma_point_A1": a1_counts[idx] // (n * point_order),
            "physical_sha256": case_metadata(shape, m)["physical_sha256"],
        }
    return result


def self_test():
    assert count_blocks((1000,))["m3"]["gamma_point_A1"] == 83333
    assert count_blocks((32, 32))["m3"]["gamma_point_A1"] == 22291
    cubic = count_blocks((10, 10, 10))
    assert cubic["m3"]["fixed_magnetization"] == 166167000
    assert cubic["m3"]["gamma_translation"] == 166167
    assert cubic["m3"]["gamma_point_A1"] == 4299
    assert cubic["m2"]["gamma_translation"] == 503
    print("PASS: independent large-site Burnside dimensions, including cubic A1g D=4299")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--shape", default="10,10,10")
    p.add_argument("--all-primary", action="store_true")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--max-sites", type=int, default=2048)
    a = p.parse_args()
    if a.self_test:
        self_test()
    else:
        shapes = PRIMARY_SHAPES if a.all_primary else (parse_shape(a.shape),)
        for shape in shapes:
            print(json.dumps(count_blocks(shape, a.max_sites), sort_keys=True))
