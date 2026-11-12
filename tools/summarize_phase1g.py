#!/usr/bin/env python3
"""Reduce Phase 1-G results without erasing failures or arithmetic-mode follow-ups."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
from numpy import __version__ as numpy_version

BASE='00702f66d404f9b81fd7db241b7b5da30fb71e2d'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(p,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,keys,lineterminator='\n');writer.writeheader()
        for row in rows:writer.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['runs','default-runs','fp-runs','output']:p.add_argument('--'+arg,type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    repo=Path(__file__).resolve().parents[1];summary={};selected={};perf_rows=[];perf_stats=[]
    manifests={key:read(root/'manifest.json') for key,root in [('required',a.runs),('default',a.default_runs),('fp_followup',a.fp_runs)]}
    for label,m in manifests.items():write(out/(label+'-commands.json'),m)
    assert all(c['exit_code']==0 for c in manifests['required']['commands'])
    assert all(c['exit_code']==0 or c['label']=='diagnostics' for c in manifests['default']['commands'])
    assert all(c['exit_code']==0 for c in manifests['fp_followup']['commands'])
    protected=['ham2vec','state_lists','input_param','blas_lapack_interfaces','main','get_expectation_values']
    protected_sha={}
    for name in protected:
        path=repo/'source'/(name+'.f90')
        assert path.read_bytes()==subprocess.check_output(['git','show',BASE+':source/'+name+'.f90'],cwd=repo)
        protected_sha[name]=sha(path)
    before=subprocess.check_output(['git','show',BASE+':source/krylov_solvers.f90'],cwd=repo,text=True)
    after=(repo/'source/krylov_solvers.f90').read_text()
    helper_names=['finite_real','normalize_vector','energy_change','residual_limit','near_limit','random_start','orthogonalize','complement']
    for name in helper_names:
        pattern=rf'^  [^\n]*\b(?:function|subroutine) {name}\(.*?^  end (?:function|subroutine) {name}'
        assert re.search(pattern,before,re.M|re.S).group()==re.search(pattern,after,re.M|re.S).group(),name
    for name in ['krylov_solvers','lanczos','eigen_solver']:
        original=subprocess.check_output(['git','show',BASE+':source/'+name+'.f90'],cwd=repo,text=True)
        current=(repo/'source'/(name+'.f90')).read_text()
        for token in ['write(*','namelist /']:
            assert [l.strip() for l in original.splitlines() if token in l]==[l.strip() for l in current.splitlines() if token in l]
    for profile in ['o2','gnu13','bounds','default']:
        root=(a.default_runs if profile=='default' else a.runs)/profile
        r=read(root/'residual-policy/report.json');assert r['before_equivalence']==146 and len(r['cases'])==146
        assert r['harness_sha256']==sha(repo/'tests/residual_policy_validation.py')
        info={k:v for k,v in r.items() if k not in ['cases','controls']}
        info['raw_report_sha256']=sha(root/'residual-policy/report.json');info['transitions']={}
        info['relative_cli_scale_checks']=0
        # Independent per-pair check with B=0 and nonzero origin, where the
        # lambda-dependent scales differ across ALG=2 pairs.
        sys.path.insert(0,str(repo/'tests'))
        from diagnostics_reproducibility import telemetry
        for alg in [1,2]:
            fields,pairs=telemetry(root/'residual-policy'/f'cli-relative-a{alg}'/'diagnostics.dat')
            for j,pair in enumerate(pairs,1):
                expected=abs(pair['ritz']+1.)
                assert abs(float(fields[f'requested_residual_scale_{j}'])-expected)<=3e-15*max(expected,1e-300)
                info['relative_cli_scale_checks']+=1
        sys.path.pop(0)
        for mode in r['summary']['policies']:
            transitions=Counter((g['legacy']['status'],g['policies'][mode]['status']) for g in r['cases'])
            info['transitions'][mode]={f'{x}->{y}':n for (x,y),n in transitions.items()}
            assert transitions[(1,0)]==0
        records=[]
        for g in r['cases']:
            legacy=g['legacy'];body=(root/'residual-policy'/(g['id']+'-legacy')/'input.dat').read_bytes().split(b'\n',1)[1]
            rng_keys=[k for k in legacy['fields'] if k.startswith('rng_initial_state_')]+['seed','random_vectors','rng_active_threads']
            for mode,result in [('legacy',legacy),*g['policies'].items()]:
                work=root/'residual-policy'/(g['id']+'-'+mode)
                assert (work/'input.dat').read_bytes().split(b'\n',1)[1]==body
                assert all(result['fields'][k]==legacy['fields'][k] for k in rng_keys)
                fields=result['fields'];pairs=result['pairs'];assessment=result.get('assessment',{})
                row={k:g[k] for k in ['id','family','alg','threads','seed','tol']}
                row.update(mode=mode,status=result['status'],reason=result['reason'],iterations=result['iterations'],
                           h_calls=result['h_calls'],seconds=result['seconds'],rss_kib=result['rss_kib'],input_sha256=result['input_sha256'],
                           legacy_observed_scale=fields['operator_scale'],legacy_energy_gate=fields['energy_gate_satisfied'],
                           inclusion=fields['lowest_state_inclusion'],
                           residuals=[x['residual'] for x in pairs],ritz=[x['ritz'] for x in pairs],
                           rayleigh=[x['rayleigh'] for x in pairs],norms=[x['norm'] for x in pairs],
                           legacy_thresholds=[x['threshold'] for x in pairs],legacy_residual_certificates=[x['certified'] for x in pairs],
                           requested_atol=fields.get('requested_residual_atol'),requested_rtol=fields.get('requested_residual_rtol'),
                           requested_B=fields.get('requested_operator_scale'),requested_origin=fields.get('requested_energy_origin'),
                           requested_scales=[fields.get(f'requested_residual_scale_{j}') for j in range(1,len(pairs)+1)],
                           requested_thresholds=[fields.get(f'requested_residual_threshold_{j}') for j in range(1,len(pairs)+1)],
                           requested_met=[fields.get(f'requested_residual_met_{j}') for j in range(1,len(pairs)+1)],
                           all_requested_met=fields.get('all_requested_residuals_met'),**assessment)
                records.append(row)
        table(out/(profile+'-cases.csv'),records)
        selected[profile]=[g for g in r['cases'] if g['id'] in ['0016-t1','0040-t2','0042-t1','0043-t1','0068-t1','0069-t1']]
        if profile=='o2':write(out/'analytic-controls.json',r['controls'])
        perf=read(root/'policy-timing/report.json')
        assert perf['harness_sha256']==sha(repo/'tools/phase1g_policy_timing.py')
        for row in perf['rows']:
            f=row['fields'];asmt=row.get('assessment',{})
            perf_rows.append(dict(profile=profile,**{k:row[k] for k in ['alg','threads','repeat','mode','status','h_calls','seconds','rss_kib']},
                                  reported_residuals=[p['residual'] for p in row['pairs']],**asmt))
        for alg in [1,2]:
            for threads in [1,2]:
                group=dict(profile=profile,alg=alg,threads=threads,modes={})
                for mode in ['legacy','loose','strict']:
                    rows=[x for x in perf['rows'] if (x['alg'],x['threads'],x['mode'])==(alg,threads,mode)]
                    group['modes'][mode]=dict(statuses=dict(Counter(str(x['status']) for x in rows)),
                        median_seconds=statistics.median(x['seconds'] for x in rows),
                        min_seconds=min(x['seconds'] for x in rows),max_seconds=max(x['seconds'] for x in rows),
                        rss_range_kib=[min(x['rss_kib'] for x in rows),max(x['rss_kib'] for x in rows)],
                        h_calls=sorted(set(x['h_calls'] for x in rows)),
                        max_reported_residual=max(max(p['residual'] for p in x['pairs']) for x in rows),
                        max_independent_residual=max(max(x.get('assessment',{}).get('residual_norms',[0.])) for x in rows))
                perf_stats.append(group)
        fp=read(a.fp_runs/profile/'results/report.json')
        assert fp['before_comparisons']==136 and fp['comparison_cases']==136
        assert fp['harness_sha256']==sha(repo/'tests/diagnostics_reproducibility.py')
        info['final_diagnostics']={k:v for k,v in fp.items() if k not in ['cases','controls','performance']}
        info['final_diagnostics']['report_sha256']=sha(a.fp_runs/profile/'results/report.json')
        info['initial_diagnostics_log_sha256']=sha(root/'diagnostics.log')
        zero=next(p for p in (a.fp_runs/profile/'results').glob('*-on') if (p/'float-mode.dat').exists())
        info['gradual_underflow']=(zero/'float-mode.dat').read_text().strip()=='T'
        suites={}
        for label in ['independent','fields','stability','stability-cli','stress','precision']:
            report=read(root/label/'report.json')
            suites[label]={k:v for k,v in report.items() if k in ['passed','failed','counts','cli_counts','operator_counts','permutation_counts']}
            suites[label]['report_sha256']=sha(root/label/'report.json')
        for label in ['standard','lapack']:
            log=(root/(label+'.log')).read_text()
            if label=='standard':assert '# TOTAL: 1' in log and '# FAIL:  0' in log
            else:assert 'PASS: 11 interface/error-path checks' in log
            suites[label]=dict(passed=1 if label=='standard' else 11,failed=0,log_sha256=sha(root/(label+'.log')))
        if profile!='default':
            h=read(root/'historical.json');assert h['passed']==27 and h['failed']==0
            write(out/(profile+'-historical.json'),h)
            assert len(re.findall(r'^PASS: test_ex[1-9]\.sh$',(root/'extended.log').read_text(),re.M))==9
            suites['extended']=dict(passed=9,failed=0,log_sha256=sha(root/'extended.log'))
            suites['historical']=dict(passed=27,failed=0,report_sha256=sha(root/'historical.json'))
            acceptance=read(root/'acceptance/report.json')
            cols=['id','alg','tol','threads','seed','exit_code','classification','eigenvalue_errors','residual_norms','max_magnetization_error','max_correlation_error']
            table(out/(profile+'-acceptance.csv'),[{k:v for k,v in row.items() if k in cols} for row in acceptance['cases']])
            failures=[x for x in acceptance['cases'] if x['classification']=='solver_failure']
            assert all(x['exit_code']!=0 and not (root/'acceptance'/x['id']/'eigenvalues.dat').exists() for x in failures)
            suites['acceptance']=dict(counts=acceptance['counts'],report_sha256=sha(root/'acceptance/report.json'),failure_eigenvalue_files=0)
        else:
            for label in ['extended','historical','acceptance']:suites[label]={'status':'unrun in supplemental quick profile'}
        info['suites']=suites;summary[profile]=info
    write(out/'summary.json',summary);write(out/'selected-cases.json',selected)
    table(out/'timing.csv',perf_rows);write(out/'timing-summary.json',perf_stats)
    libraries={}
    for manifest in [manifests['required'],manifests['default']]:
        for profile,spec in manifest['profiles'].items():
            for name,digest in spec['after_source_sha256'].items():assert sha(repo/'source'/name)==digest
            for lib in re.findall(r'=> (/[^ ]+)',spec['linked_libraries']):libraries[str(Path(lib).resolve())]=sha(Path(lib).resolve())
    write(out/'provenance.json',dict(base=BASE,protected_source_sha256=protected_sha,unchanged_legacy_helpers=helper_names,
        unchanged_normal_stdout_and_namelists=True,linked_libraries_sha256=libraries,
        observed_after_validation=dict(cpu=next(l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),
             platform=platform.platform(),python=sys.version,numpy=numpy_version,affinity=sorted(os.sched_getaffinity(0)),
             cpu_max=Path('/sys/fs/cgroup/cpu.max').read_text().strip()),
        source_sha256={str(p.relative_to(repo)):sha(p) for p in [*sorted((repo/'source').glob('*.f90')),
             repo/'tests/krylov_policy_driver.F90',repo/'tests/residual_policy_validation.py',repo/'tests/krylov_diagnostics_driver.F90',
             repo/'tests/diagnostics_reproducibility.py',repo/'tools/phase1g_validation.py',repo/'tools/phase1g_policy_timing.py',
             repo/'tools/phase1g_fp_validation.py',Path(__file__)]}))
    print(json.dumps({k:v['summary'] for k,v in summary.items()}))


if __name__=='__main__':main()
