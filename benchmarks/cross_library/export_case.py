#!/usr/bin/env python3
"""Export canonical U(1)-conserving NN Heisenberg bond list and metadata."""
import argparse
import hashlib
import json
from pathlib import Path
from reference import bonds_for_shape, case_metadata, parse_shape


def export_case(shape, ndown, output_prefix):
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    meta = case_metadata(shape, ndown)
    if meta["sites"] > 100000:
        raise ValueError("site count exceeds deliberate export safety bound")
    bonds = bonds_for_shape(shape)
    text = "".join("%d %d\n" % (i, j) for i, j in bonds)
    bpath = prefix.with_suffix(".bonds")
    mpath = prefix.with_suffix(".json")
    bpath.write_text(text, encoding="ascii")
    manifest = dict(meta)
    manifest.update(
        bond_file=bpath.name, bond_sha256=hashlib.sha256(text.encode("ascii")).hexdigest(),
        bond_convention="0-based sorted undirected NN, one line 'i j', no header",
        model_terms={"zz": 1.0, "+-": 0.5, "-+": 0.5},
        safety="data only; no Hilbert-space allocation or eigensolver invoked")
    mpath.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shape", required=True)
    p.add_argument("--ndown", type=int, required=True)
    p.add_argument("--output-prefix", required=True)
    args = p.parse_args()
    print(json.dumps(export_case(parse_shape(args.shape), args.ndown,
                                 args.output_prefix), sort_keys=True))


if __name__ == "__main__":
    main()
