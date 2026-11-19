#!/usr/bin/env python3
"""Audit native L4/L6 current-main DM table against frozen CI-produced report.

This is an editorial transcription/provenance check, NOT a re-run of ED2
nor direct residual certification. The end-to-end Fortran verification is
a separate CI job. Pure Python stdlib.
"""
import argparse
import json
import math
from pathlib import Path

START="% MAINLINE_CURRENT_NATIVE_DM_TABLE_BEGIN:"
END="% MAINLINE_CURRENT_NATIVE_DM_TABLE_END"


def verify(report,tex):
    if report.get("status")!="review_only_native_reproduction_not_release":
        raise ValueError("Missing explicit provisional status")
    if report.get("no_direct_residual_sidecar") is not True:
        raise ValueError("Original mainline residual caveat not encoded")
    if report.get("native_artifact_id")!=14492907330:
        raise ValueError("Unexpected native provenance artifact")
    if "tab:mainline-current-native-dm" not in tex:
        raise ValueError("Current-main DM table missing")
    if tex.count(START)!=1 or tex.count(END)!=1:
        raise ValueError("Expected exactly one source-audited data table")
    block=tex.split(START,1)[1].split(END,1)[0]
    cases=report.get("sizes")
    if not isinstance(cases,list) or len(cases)!=2:
        raise ValueError("Expected two native checked system sizes")
    for case,n in zip(cases,(4,6)):
        if case["sites"]!=n or case["dimension"]!=6**(n//2):
            raise ValueError("Incorrect physical full-space size")
        if case["run_count"]!=3 or case["direct_residual_certified"] is not False:
            raise ValueError("Missing required three-run/native sidecar scope")
        fields=("energy_e0_dm0","chirality_k_dm0",
                "hf_abs_error","chi_per_bond",
                "independent_reference_chi_per_bond")
        if not all(math.isfinite(float(case[k])) for k in fields):
            raise ValueError("Non-finite scientific data")
        if case["hf_abs_error"]>3e-5:
            raise ValueError("Hellmann--Feynman identity mismatch")
        if abs(case["chi_per_bond"]-
               case["independent_reference_chi_per_bond"])>3e-3:
            raise ValueError("Susceptibility reference mismatch")
        expected=(f'{n} & {case["energy_e0_dm0"]:.12f} & '
                  f'{case["chirality_k_dm0"]:.12f} & '
                  f'{case["hf_abs_error"]:.2e} & '
                  f'{case["chi_per_bond"]:.9f}')
        if expected not in block or block.count(expected)!=1:
            raise ValueError("Native response table/JSON mismatch: "+expected)
    return cases


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,
                        default=Path("paper/data/current_main_native_dm_l4_l6_20261119.json"))
    parser.add_argument("--tex",type=Path,default=Path("paper/cpc_ed2_main.tex"))
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    data=json.loads(args.json.read_text(encoding="utf-8"))
    tex=args.tex.read_text(encoding="utf-8")
    verify(data,tex)
    if args.self_test:
        original='0.568732921'
        if tex.count(original)!=1:
            raise ValueError("Missing unique six-site native susceptibility")
        try:
            verify(data,tex.replace(original,'0.668732921'))
        except ValueError:
            print("CPC_NATIVE_DM_TABLE_NEGATIVE_TEST_PASS")
        else:
            raise AssertionError("Corrupted 6-site susceptibility was not detected")
    print("CPC_CURRENT_MAIN_DM_TABLE_PASS 2/2 native source-linked L4/L6 rows")


if __name__=="__main__":
    main()
