#!/usr/bin/env python3
"""Non-executing large-sector capacity and symmetry preflight.

No ED2 diagonalization, no PBS job, and no assumption that a 1000-site
fixed-magnon calculation fits a given memory allocation.
"""
import argparse
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import cpc_pointgroup_profile as pg

def estimate(side):
    assert side % 2 == 0 and side >= 4
    n = side * side
    group = 8 * n
    images = group * (n + 1)
    raw = math.comb(n, 3)
    result = dict(side=side, sites=n, raw_three_magnon_dimension=raw,
                  c4v_group_elements=group, site_image_entries=images,
                  permutation_bytes_if_int32=4 * images,
                  group_count_limit=8192, image_entry_limit=12000000,
                  eligible_by_validation_caps=(group <= 8192 and images <= 12000000))
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", default="12,16,20,26,32")
    parser.add_argument("--independent-projector", action="store_true",
                        help="Compute character/Burnside dimensions; can be slow at 1024 sites")
    args = parser.parse_args()
    import json
    for side in map(int, args.sizes.split(",")):
        row = estimate(side)
        if args.independent_projector:
            row["independent_dimensions"] = {
                sector: pg.independent_dimension(side*side, 3, pg.symmetry_generators(side, sector))
                for sector in ("gamma-t", "gamma-a1", "gamma-b1", "x-a1")}
        print(json.dumps(row, sort_keys=True))
    assert estimate(32)["eligible_by_validation_caps"]
    assert not estimate(34)["eligible_by_validation_caps"]
if __name__ == "__main__":
    main()
