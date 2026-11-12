#!/usr/bin/env python3
"""Reduce final Phase 1-F runs; retain initial failures and validated corrections.

Requires the full runner's output plus fresh stability-format-fixed{,-cli} runs
for the test-only exponent-format correction. Never overwrites an evidence dir.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
import platform
import sys
from pathlib import Path
import re
import statistics
import subprocess
from numpy import __version__ as numpy_version

BASE='b756b0009f08027e8a2ae8a1e349d78d8e81d9c5'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def write(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def table(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fields,lineterminator='\n');w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in row.items()})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--success-timing',type=Path,required=True);a=p.parse_args()
    root=a.runs.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    repo=Path(__file__).resolve().parents[1];manifest=read(root/'manifest.json');write(out/'command-manifest.json',manifest)
    old_driver=subprocess.check_output(['git','show',BASE+':tests/lanczos_stability_driver.f90'],cwd=repo,text=True)
    corrected=(repo/'tests/lanczos_stability_driver.f90').read_text()
    assert old_driver.replace('es25.15','es26.15e3')==corrected
    protected=['ham2vec','state_lists','input_param','blas_lapack_interfaces','main','get_expectation_values']
    protected_sha={}
    for name in protected:
        path=repo/'source'/(name+'.f90')
        original=subprocess.check_output(['git','show',BASE+':source/'+name+'.f90'],cwd=repo)
        assert path.read_bytes()==original,name
        protected_sha[name]=sha(path)
    old_kernel=subprocess.check_output(['git','show',BASE+':source/krylov_solvers.f90'],cwd=repo,text=True)
    new_kernel=(repo/'source/krylov_solvers.f90').read_text()
    functions=['finite_real','normalize_vector','energy_change','residual_limit','near_limit','orthogonalize','complement']
    for name in functions:
        pattern=rf'^  [^\n]*\b(?:function|subroutine) {name}\(.*?^  end (?:function|subroutine) {name}'
        old=re.search(pattern,old_kernel,re.M|re.S).group();new=re.search(pattern,new_kernel,re.M|re.S).group()
        assert old==new,name
    # Normal production stdout/namelist statements have no added fields/formats.
    for name in ['krylov_solvers','lanczos','eigen_solver']:
        old=subprocess.check_output(['git','show',BASE+':source/'+name+'.f90'],cwd=repo,text=True)
        new=(repo/'source'/(name+'.f90')).read_text()
        for marker in ['write(*','namelist /']:
            assert [line.strip() for line in old.splitlines() if marker in line]==[line.strip() for line in new.splitlines() if marker in line]
    results={};performance=[];perf_summary=[];comparison=[];selected={};resolved=[]
    for profile in ['o2','gnu13','bounds']:
        pr=root/profile;d=read(pr/'diagnostics/report.json')
        assert d['comparison_cases']==d['before_comparisons']==136
        assert len(d['seeded_replays'])==24 and len(d['cli_checks'])==22 and len(d['error_paths'])==6
        assert d['harness_sha256']==sha(repo/'tests/diagnostics_reproducibility.py')
        for name,digest in manifest['profiles'][profile]['after_source_sha256'].items():
            assert sha(repo/'source'/name)==digest,(profile,name)
        info={k:v for k,v in d.items() if k not in ['cases','performance','controls']}
        info['raw_diagnostics_report_sha256']=sha(pr/'diagnostics/report.json')
        for case in d['cases']:
            fields=case['diagnostics']['fields']
            comparison.append(dict(profile=profile,id=case['id'],family=case['family'],status=case['status'],
                                   reason=case['reason'],h_calls=case['h_calls'],vector_sha256=case['vector_sha256'],
                                   fields=fields,pairs=case['diagnostics']['pairs']))
        selected[profile]=[case for case in d['cases'] if case['id'] in ['0017-t1','0040-t2','0042-t1','0043-t1']]
        if profile=='o2':write(out/'analytic-controls.json',d['controls'])
        for row in d['performance']:performance.append(dict(profile=profile,**row))
        for alg in [1,2]:
            modes={m:[r for r in d['performance'] if r['alg']==alg and r['mode']==m] for m in ['before','off','on']}
            stats={m:dict(solver_statuses=sorted(set(r['status'] for r in rows)),median_seconds=statistics.median(r['seconds'] for r in rows),
                          min_seconds=min(r['seconds'] for r in rows),max_seconds=max(r['seconds'] for r in rows),
                          min_rss_kib=min(r['rss_kib'] for r in rows),max_rss_kib=max(r['rss_kib'] for r in rows),
                          h_calls=sorted(set(r['h_calls'] for r in rows))) for m,rows in modes.items()}
            assert stats['before']['h_calls']==stats['off']['h_calls']==stats['on']['h_calls']
            perf_summary.append(dict(profile=profile,alg=alg,workload='original-cap12',**stats,
                                     on_over_off=stats['on']['median_seconds']/stats['off']['median_seconds'],
                                     off_over_before=stats['off']['median_seconds']/stats['before']['median_seconds']))
        suites={}
        for label in ['independent','fields','stress','precision']:
            report=read(pr/label/'report.json')
            suites[label]={k:v for k,v in report.items() if k in ['passed','failed','counts','cli_counts','operator_counts','permutation_counts']}
            suites[label]['report_sha256']=sha(pr/label/'report.json')
        for label in ['stability','stability-cli']:
            if not (pr/label/'report.json').exists():
                suites['initial-'+label]={'status':'unrun after preceding suite failure'}
                continue
            report=read(pr/label/'report.json')
            suites['initial-'+label]=dict(passed=report['passed'],failed=report['failed'],report_sha256=sha(pr/label/'report.json'))
            if report['failed']:
                resolved.append(dict(profile=profile,suite=label,
                                     original_failures=[case for case in report['cases'] if case['status']=='failed']))
        for label in ['stability-format-fixed','stability-format-fixed-cli']:
            report=read(pr/label/'report.json');assert report['failed']==0
            suites[label]=dict(passed=report['passed'],failed=report['failed'],report_sha256=sha(pr/label/'report.json'))
            write(out/(profile+'-'+label+'.json'),report)
        historical=read(pr/'historical.json');assert historical['passed']==27 and historical['failed']==0
        write(out/(profile+'-historical.json'),historical)
        suites['historical']=dict(passed=27,failed=0,report_sha256=sha(pr/'historical.json'))
        extended=(pr/'extended.log').read_text()
        assert len(re.findall(r'^PASS: test_ex[1-9]\.sh$',extended,re.M))==9
        suites['extended']=dict(passed=9,failed=0,log_sha256=sha(pr/'extended.log'))
        assert '# TOTAL: 1' in (pr/'standard.log').read_text() and '# FAIL:  0' in (pr/'standard.log').read_text()
        assert 'PASS: 11 interface/error-path checks' in (pr/'lapack.log').read_text()
        suites['standard']=dict(passed=1,failed=0,log_sha256=sha(pr/'standard.log'))
        suites['lapack']=dict(passed=11,failed=0,log_sha256=sha(pr/'lapack.log'))
        acceptance=read(pr/'acceptance/report.json');assert len(acceptance['cases'])==1024
        failed=[case for case in acceptance['cases'] if case['classification']=='solver_failure']
        assert all(case['exit_code']!=0 and not (pr/'acceptance'/case['id']/'eigenvalues.dat').exists() for case in failed)
        assert all('maximum restart iterations without convergence' in (pr/'acceptance'/case['id']/'stdout.log').read_text() for case in failed)
        fields=['id','alg','tol','threads','seed','exit_code','classification','eigenvalue_errors','residual_norms',
                'max_magnetization_error','max_correlation_error','orthogonality_error']
        table(out/(profile+'-acceptance.csv'),[{k:v for k,v in row.items() if k in fields} for row in acceptance['cases']])
        suites['acceptance']=dict(counts=acceptance['counts'],report_sha256=sha(pr/'acceptance/report.json'),
                                  failures_are_budget_exhaustion=True,failure_eigenvalue_files=0,
                                  executable_sha256=acceptance['executable_sha256'],harness_sha256=acceptance['harness_sha256'])
        probe=pr/'probe'
        write(out/(profile+'-probe-provenance.json'),dict(provenance=read(probe/'provenance.json'),
                                                        instrumentation_lines=(probe/'phase1c-instrumentation.patch').read_text().splitlines()))
        info['suites']=suites;results[profile]=info
    # The raw manifest deliberately keeps original parse failures. Only that
    # diagnosed suite has a validated replacement report; no solver failure is
    # turned into success here.
    for cmd in manifest['commands']:
        assert cmd['exit_code']==0 or cmd['label']=='stability',cmd
    supplement=read(a.success_timing/'report.json')
    assert supplement['harness_sha256']==sha(repo/'tools/phase1f_success_timing.py')
    assert all(row['status']==0 for row in supplement['rows'])
    write(out/'converged-alg2-timing.json',supplement)
    for profile in ['o2','gnu13','bounds']:
        modes={m:[r for r in supplement['rows'] if r['profile']==profile and r['mode']==m] for m in ['before','off','on']}
        stats={m:dict(solver_statuses=[0],median_seconds=statistics.median(r['seconds'] for r in rows),
                      min_seconds=min(r['seconds'] for r in rows),max_seconds=max(r['seconds'] for r in rows),
                      min_rss_kib=min(r['rss_kib'] for r in rows),max_rss_kib=max(r['rss_kib'] for r in rows),
                      h_calls=sorted(set(r['h_calls'] for r in rows))) for m,rows in modes.items()}
        assert stats['before']['h_calls']==stats['off']['h_calls']==stats['on']['h_calls']
        perf_summary.append(dict(profile=profile,alg=2,workload='converged-cap32',**stats,
                                 on_over_off=stats['on']['median_seconds']/stats['off']['median_seconds'],
                                 off_over_before=stats['off']['median_seconds']/stats['before']['median_seconds']))
    write(out/'summary.json',results);write(out/'performance-summary.json',perf_summary)
    table(out/'performance.csv',performance);table(out/'kernel-comparisons.csv',comparison)
    write(out/'selected-observations.json',selected)
    write(out/'format-correction.json',dict(initial_failures=resolved,
                                           change='test-only ES25.15 to ES26.15E3; scientific assertions unchanged',
                                           corrected_source_sha256=sha(repo/'tests/lanczos_stability_driver.f90')))
    libraries={}
    for profile in manifest['profiles'].values():
        for path in re.findall(r'=> (/[^ ]+)',profile['linked_libraries']):
            library=Path(path).resolve()
            libraries[str(library)]=sha(library)
    machine=dict(observed_after_validation=True,platform=platform.platform(),python=sys.version,numpy=numpy_version,
                 cpu=next(line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                          if line.startswith('model name')),affinity=sorted(os.sched_getaffinity(0)),
                 cpu_max=Path('/sys/fs/cgroup/cpu.max').read_text().strip(),
                 thread_environment=dict(OMP_NUM_THREADS=[1,2],OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'))
    write(out/'provenance.json',dict(base=BASE,machine=machine,linked_library_sha256=libraries,
                                    protected_source_sha256=protected_sha,unchanged_numerical_helpers=functions,
                                    production_stdout_and_namelists_unchanged=True,
                                    maintained_source_sha256={str(path.relative_to(repo)):sha(path) for path in
                                    [*sorted((repo/'source').glob('*.f90')),repo/'tests/krylov_diagnostics_driver.F90',
                                     repo/'tests/diagnostics_reproducibility.py',repo/'tools/phase1f_validation.py',
                                     repo/'tools/phase1f_success_timing.py',Path(__file__)]}))
    print(json.dumps({p:dict(before_comparisons=r['before_comparisons'],statuses=r['diagnostic_solver_statuses'],
                             historical=r['suites']['historical']['passed'],acceptance=r['suites']['acceptance']['counts'])
                      for p,r in results.items()}))


if __name__=='__main__':main()
