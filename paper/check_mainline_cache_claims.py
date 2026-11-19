#!/usr/bin/env python3
"""Check CPC mainline generic cache claims against source and exact arithmetic.

This tests the numerical *transcription* and published flag provenance.
The model-dependent Nrec=191808 was measured in prior native local
mixed-spin real/complex DM calculations; it is NOT independently derived
from this Python program. Do not call this a runtime/RSS benchmark.
"""
import argparse
from pathlib import Path
import re

D = 7776
PADDED_WIDTH = 93
CERTIFIED_WIDTH = 57
NREC = 191808
BATCH_ROWS = 64
ENTRY_BYTES = 4 + 16


def expected():
    padded = ENTRY_BYTES * PADDED_WIDTH * D
    bounded = ENTRY_BYTES * CERTIFIED_WIDTH * D
    csr = 4 * (D + 1) + ENTRY_BYTES * NREC
    stage = ENTRY_BYTES * CERTIFIED_WIDTH * BATCH_ROWS
    return {
        "Original unbounded padded representation": padded,
        "Coefficient-certified padded representation": bounded,
        "Coefficient-bounded padded + CSR simultaneously": bounded + csr,
        "Two-pass direct CSR + staging (upper bound)": csr + stage,
    }


def verify(tex, source):
    for name in (
        "ED2_CACHE_COEFF_BOUND", "ED2_CSR_COMPLEX_CACHE",
        "ED2_CSR_RELEASE_PADDED_CACHE", "ED2_CSR_DIRECT_BUILD",
    ):
        if source.count("'" + name + "'") < 1:
            raise ValueError("Required default-off mainline opt-in missing: " + name)
    if "integer, parameter :: direct_batch_rows=64" not in source:
        raise ValueError("Manuscript batch width no longer matches generic main")
    if r"\label{tab:general-csr-payload}" not in tex:
        raise ValueError("Missing general cache-memory table")
    if r"\appendix" not in tex or r"\section{Unmerged cubic-specific" not in tex:
        raise ValueError("Unmerged lattice-specific trial must be appendix")
    for title, bytes_expected in expected().items():
        pat = (re.escape(title) + r"\s*&\s*"
               + format(bytes_expected, ",").replace(",", r",")
               + r"\s*\\{2}")
        if not re.search(pat, tex):
            raise ValueError("General cache memory table differs: "
                             + title + ", expected " + format(bytes_expected, ","))
    if "not total process" not in tex and "not whole-process" not in tex:
        raise ValueError("Missing allocation vs RSS caveat")
    return expected()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tex", type=Path, default=Path("paper/cpc_ed2_main.tex"))
    p.add_argument("--source", type=Path, default=Path("source/ham2vec.f90"))
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    tex = args.tex.read_text(encoding="utf-8")
    src = args.source.read_text(encoding="utf-8")
    results = verify(tex, src)
    if args.self_test:
        before = format(results["Two-pass direct CSR + staging (upper bound)"], ",")
        try:
            verify(tex.replace(before, "9,999,999", 1), src)
        except ValueError:
            pass
        else:
            raise AssertionError("Checker accepted corrupted memory figure")
        print("CPC_GENERIC_CACHE_CLAIM_NEGATIVE_TEST_PASS")
    print("CPC_GENERIC_CACHE_CLAIM_PASS D={} Nrec={} Bstage={} table={}"
          .format(D, NREC, BATCH_ROWS, results))


if __name__ == "__main__":
    main()
