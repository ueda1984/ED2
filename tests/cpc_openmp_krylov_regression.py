#!/usr/bin/env python3
"""Compare numerical outcomes without assuming identical Krylov trajectories.

No solver tolerance or independent accuracy threshold is relaxed. All raw
results, warnings, nonconvergences and iteration differences are retained.
An early audit incorrectly required iteration counts to match across different
compiler flags/BLAS GEMM widths. Accurate eigenpairs can require a different
number of restarts. This gate instead requires the existing independent
energy/subspace/residual/orthogonality assessments and solver classifications
to agree. Trajectory changes remain visible, not silently discarded.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harness',required=True,type=Path)
    parser.add_argument('--baseline-driver',required=True,type=Path)
    parser.add_argument('--candidate-driver',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--threads',default='1,2,4')
    args=parser.parse_args()
    threads=tuple(int(x) for x in args.threads.split(','))
    if not threads or min(threads)<1 or len(set(threads))!=len(threads):
        parser.error('Invalid thread counts')
    sys.path.insert(0,str(args.harness.resolve()))
    from krylov_stress import definitions,run
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    for engine in ('baseline','candidate'):(root/engine).mkdir()
    old=args.baseline_driver.resolve(strict=True)
    new=args.candidate_driver.resolve(strict=True)
    fixtures=list(definitions());records=[];failures=[];differences=[]
    categorical=('classification','solver_status','exit_code',
                 'strict_energy_target_met','requested_space_target_met',
                 'valid_pair_under_loose_assessment')
    trajectory=('iterations','restart_cycles','h_applications')
    def save():
        result=dict(passed=len(records)==len(fixtures)*len(threads) and not failures,
                    fixture_count=len(fixtures),thread_counts=threads,
                    expected_comparisons=len(fixtures)*len(threads),completed=len(records),
                    failures=failures,trajectory_changes=differences,
                    baseline_counts=dict(Counter(x['baseline']['classification'] for x in records)),
                    candidate_counts=dict(Counter(x['candidate']['classification'] for x in records)),
                    baseline_binary_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
                    candidate_binary_sha256=hashlib.sha256(new.read_bytes()).hexdigest(),
                    records=records,
                    note='No accuracy thresholds changed; iteration counts are measured rather than required equal.')
        (root/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        compact={k:v for k,v in result.items() if k!='records'}
        (root/'summary.json').write_text(json.dumps(compact,indent=2,allow_nan=False)+'\n')
        return result
    for index,fixture in enumerate(fixtures):
        for count in threads:
            a=run(old,root/'baseline',fixture,count,index)
            b=run(new,root/'candidate',fixture,count,index)
            errors={k:[a.get(k),b.get(k)] for k in categorical if a.get(k)!=b.get(k)}
            for key in ('input_sha256','matrix_sha256','rng_state'):
                if a.get(key)!=b.get(key):errors[key]='provenance mismatch'
            if 'harness_or_runtime_error' in (a['classification'],b['classification']):
                errors['runtime']='Harness or numerical executable failed unexpectedly'
            # Explicitly preserve the independent assessment flags, including
            # the known adversarial inclusion-warning and nonconvergent cases.
            if a.get('solver_status')==0 and b.get('solver_status')==0:
                if not b.get('valid_pair_under_loose_assessment',False):
                    errors['invalid_pair']='Independent residual/orthogonality check failed'
            changed={k:[a.get(k),b.get(k)] for k in trajectory if a.get(k)!=b.get(k)}
            row=dict(index=index,threads=count,baseline=a,candidate=b,errors=errors,trajectory_changes=changed)
            records.append(row)
            if errors:failures.append(dict(index=index,threads=count,errors=errors))
            if changed:differences.append(dict(index=index,threads=count,changes=changed,
                                               baseline_classification=a['classification'],
                                               candidate_classification=b['classification']))
        if index%20==0:
            save();print('KRYLOV GATE',index,'fixtures,',len(failures),'violations,',len(differences),'trajectory changes',flush=True)
    report=save()
    print('KRYLOV GATE COMPLETE',report['completed'],'/',report['expected_comparisons'],
          'violations',len(failures),'trajectory changes',len(differences),'passed',report['passed'],flush=True)
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
