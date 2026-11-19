#!/usr/bin/env python3
"""Conservative capacity estimate for XDiag's explicit PermutationGroup API.

Pinned source: awietek/xdiag @ 27065537a743c642c35efa0e1b5ae4f326fd984c
  xdiag/symmetries/permutation_group.cpp
The implementation forms arma::Mat<int64_t> multiply(G,G), a matrix of all
G site permutations (N,G), and an inverse index vector G. This tool reports
those ALLOCATED PAYLOAD lower bounds only, never a measured peak RSS.

For cubic Gamma/Oh A1g, G = N*48 (translation semidirect Oh). The ordinary
RepresentationSet permits only one SitePermutation group; it cannot take
separate noncommuting point+translation generators as independent keys.

It is NOT a conclusion that all XDiag variants/rewrites are impossible:
a dedicated factorized-group algorithm could avoid this memory footprint.
"""
import argparse
import json
import math


def memory_floor(nsites, point_order=48, scalar_bytes=8):
    if nsites < 1 or point_order < 1:
        raise ValueError("positive site and group counts required")
    g = nsites * point_order
    multiply = scalar_bytes * g * g
    permutations = scalar_bytes * g * nsites
    inverse = scalar_bytes * g
    total = multiply + permutations + inverse
    return {
        "sites": nsites,
        "point_group_order": point_order,
        "space_group_order": g,
        "multiplication_table_bytes": multiply,
        "permutation_matrix_bytes": permutations,
        "inverse_index_bytes": inverse,
        "minimum_structural_payload_bytes": total,
        "minimum_structural_payload_gib": round(total/(1024**3), 6),
        "method": "Static source-level allocation accounting; no XDiag runtime allocation",
        "scope": "Explicit full PermutationGroup API, not other possible representations",
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--sites", type=int, default=1000)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        x = memory_floor(1000)
        assert x["space_group_order"] == 48000
        assert x["multiplication_table_bytes"] == 18432000000
        assert x["permutation_matrix_bytes"] == 384000000
        assert x["inverse_index_bytes"] == 384000
        assert x["minimum_structural_payload_bytes"] == 18816384000
        print("PASS: source-level group-table allocation accounting")
    else:
        print(json.dumps(memory_floor(a.sites), sort_keys=True))
