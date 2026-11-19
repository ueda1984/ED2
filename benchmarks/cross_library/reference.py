#!/usr/bin/env python3
"""Small dense U(1)-sector oracle for periodic spin-half Heisenberg NN benchmarks."""
import argparse
import hashlib
import itertools
import json
import math


def parse_shape(value):
    dims = tuple(int(x) for x in value.split(","))
    if not 1 <= len(dims) <= 3 or any(x < 3 for x in dims):
        raise ValueError("periodic axes must have length >=3")
    return dims


def bonds_for_shape(shape):
    strides = [math.prod(shape[:i]) for i in range(len(shape))]
    result = set()
    for i in range(math.prod(shape)):
        for axis, side in enumerate(shape):
            stride = strides[axis]
            j = i + stride if (i // stride) % side < side-1 else i - (side-1)*stride
            result.add(tuple(sorted((i, j))))
    assert len(result) == math.prod(shape)*len(shape)
    return sorted(result)


def case_metadata(shape, ndown):
    n = math.prod(shape)
    if not 0 <= ndown <= n:
        raise ValueError("invalid down-spin count")
    physical = dict(model="periodic_spin_half_Heisenberg_NN", J=1.0,
                    field=0.0, shape=list(shape), ndown=ndown,
                    site_order="x-fastest", bonds=bonds_for_shape(shape))
    digest = hashlib.sha256(json.dumps(physical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return dict(shape=list(shape), sites=n, ndown=ndown, dimension_u1=math.comb(n, ndown),
                bonds=len(physical["bonds"]), physical_sha256=digest,
                model=physical["model"], J=1.0, field=0.0, site_order="x-fastest")


def dense_reference(shape, ndown, max_dim=1024):
    import numpy as np
    meta = case_metadata(shape, ndown)
    if meta["dimension_u1"] > max_dim:
        raise ValueError("dense sector exceeds resource limit")
    masks = [sum(1 << i for i in spins) for spins in itertools.combinations(range(meta["sites"]), ndown)]
    index = {mask: j for j, mask in enumerate(masks)}
    h = np.zeros((len(masks), len(masks)), dtype=float)
    for col, mask in enumerate(masks):
        for i, j in bonds_for_shape(shape):
            opposite = ((mask >> i) ^ (mask >> j)) & 1
            h[col, col] += -0.25 if opposite else 0.25
            if opposite:
                h[index[mask ^ (1 << i) ^ (1 << j)], col] += 0.5
    assert np.array_equal(h, h.T)
    values, vectors = np.linalg.eigh(h)
    residual = np.linalg.norm(h @ vectors[:, 0] - values[0] * vectors[:, 0])
    return dict(meta, dimension=len(masks), E0=float(values[0]), direct_residual=float(residual))


def self_test():
    assert len(bonds_for_shape((6,))) == 6
    assert len(bonds_for_shape((4, 4))) == 32
    assert len(bonds_for_shape((4, 4, 4))) == 192
    assert len(bonds_for_shape((10, 10, 10))) == 3000
    assert case_metadata((10, 10, 10), 3)["dimension_u1"] == 166167000
    assert abs(dense_reference((6,), 0)["E0"] - 1.5) < 1e-12
    assert abs(dense_reference((6,), 1)["E0"] + 0.5) < 1e-12
    assert abs(dense_reference((4, 4), 1)["E0"] - 4.0) < 1e-12
    assert dense_reference((6,), 2)["direct_residual"] < 1e-12
    print("PASS: bond multiplicities, dimensions, energy and residual")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", default="6")
    ap.add_argument("--ndown", type=int, default=2)
    ap.add_argument("--max-dim", type=int, default=1024)
    ap.add_argument("--metadata-only", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
    else:
        shape = parse_shape(a.shape)
        obj = case_metadata(shape, a.ndown) if a.metadata_only else dense_reference(shape, a.ndown, a.max_dim)
        print(json.dumps(obj, sort_keys=True))
