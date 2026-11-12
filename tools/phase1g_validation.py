#!/usr/bin/env python3
"""Reproduce Phase 1-G validation in fresh external trees; never build the checkout.

Cloud example: activate /workspace/ed2-env/activate.sh, then run with --output NEW.
GNU 13's path is configurable. All commands/statuses and snapshot hashes are kept;
a failed test is recorded and independent checks continue. No worktrees are used.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

BASE='00702f66d404f9b81fd7db241b7b5da30fb71e2d'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--profiles',nargs='+',choices=['o2','gnu13','bounds','default'],default=['o2','gnu13','bounds'])
    p.add_argument('--gnu13',default='/workspace/ed2-env/gnu13/usr/bin/gfortran-13')
    p.add_argument('--ldflags',default='-L/workspace/ed2-env/tools/usr/lib/x86_64-linux-gnu/openblas-pthread -L/workspace/ed2-env/tools/usr/lib/gcc/x86_64-linux-gnu/14')
    p.add_argument('--quick',action='store_true',help='Omit extended/historical/seeded-probe suites; record omissions')
    a=p.parse_args();repo=a.repo.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source_files=subprocess.check_output(['git','ls-files','-z','--cached','--others','--exclude-standard'],cwd=repo).decode().split('\0')
    after_archive=out/'after-source.tar'
    with tarfile.open(after_archive,'w') as archive:
        for name in sorted(set(filter(None,source_files))):
            if (repo/name).exists():archive.add(repo/name,arcname=name,recursive=False)
    before_archive=out/'before-source.tar'
    before_archive.write_bytes(subprocess.check_output(['git','archive',BASE],cwd=repo))
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for key in ['ED2_RNG_SEED','ED2_INITIAL_VECTOR','ED2_DIAGNOSTICS','ED2_RESIDUAL_ATOL','ED2_RESIDUAL_RTOL','ED2_RESIDUAL_SCALE','ED2_ENERGY_ORIGIN']:env.pop(key,None)
    manifest=dict(base=BASE,head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                  working_tree=subprocess.check_output(['git','status','--short'],cwd=repo,text=True),
                  archives={p.name:sha(p) for p in [before_archive,after_archive]},
                  harness_sha256=sha(Path(__file__)),quick=a.quick,commands=[],profiles={})
    def save():
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def run(profile,label,cmd,cwd,required=True,timeout=2400):
        log=out/profile/(label+'.log');start=time.perf_counter()
        with log.open('w') as stream:
            try:
                r=subprocess.run(list(map(str,cmd)),cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
                code=r.returncode
            except subprocess.TimeoutExpired:
                code=124
        manifest['commands'].append(dict(profile=profile,label=label,args=list(map(str,cmd)),cwd=str(cwd),
                                         exit_code=code,seconds=time.perf_counter()-start,required=required,log_sha256=sha(log)))
        save();print(profile,label,'exit',code,flush=True)
        return code==0
    save()
    for profile in a.profiles:
        root=out/profile;root.mkdir()
        flags='-O2 -fopenmp' if profile!='bounds' else '-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra'
        if profile=='default':flags='-fopenmp -Ofast'
        fc=a.gnu13 if profile=='gnu13' else 'gfortran'
        builds={}
        for label,archive_path in [('before',before_archive),('after',after_archive)]:
            stage=root/label;stage.mkdir()
            with tarfile.open(archive_path) as archive:archive.extractall(stage,filter='data')
            okay=run(profile,label+'-autoreconf',['autoreconf','-vfi'],stage)
            okay=okay and run(profile,label+'-configure',['./configure','FC='+fc,'FCFLAGS='+flags,'LDFLAGS='+a.ldflags],stage)
            okay=okay and run(profile,label+'-clean',['make','clean'],stage)
            okay=okay and run(profile,label+'-build',['make','-j2'],stage)
            if not okay:break
            builds[label]=stage
        if len(builds)!=2:continue
        before=builds['before'];stage=builds['after']
        manifest['profiles'][profile]=dict(compiler=subprocess.check_output([fc,'--version'],text=True).splitlines()[0],
                                           flags=flags,ldflags=a.ldflags,
                                           after_source_sha256={p.name:sha(p) for p in (stage/'source').glob('*.f90')},
                                           before_source_sha256={p.name:sha(p) for p in (before/'source').glob('*.f90')},
                                           linked_libraries=subprocess.check_output(['ldd',str(stage/'source/QS3ED2')],env=env,text=True))
        run(profile,'diagnostics-driver',['make','-C','tests','krylov_diagnostics_driver','krylov_policy_driver'],stage)
        before_driver=root/'before-diagnostics-driver'
        before_policy=root/'before-policy-driver'
        for label,source,binary in [('before-driver','krylov_diagnostics_driver.F90',before_driver),
                                    ('before-policy-driver','krylov_policy_driver.F90',before_policy)]:
            run(profile,label,[fc,*flags.split(),'-cpp','-DPHASE1F','-I'+str(before/'source'),stage/'tests'/source,
                               before/'source/krylov_solvers.o',before/'source/krylov_observations.o',
                               before/'source/blas_lapack_interfaces.o',*a.ldflags.split(),
                               '-llapack','-lopenblas','-o',binary],root)
        run(profile,'residual-policy',[sys.executable,'-B',stage/'tests/residual_policy_validation.py',
                                      '--driver',stage/'tests/krylov_policy_driver','--before-driver',before_policy,
                                      '--exe',stage/'source/QS3ED2','--output',root/'residual-policy'],stage)
        run(profile,'policy-timing',[sys.executable,'-B',stage/'tools/phase1g_policy_timing.py',
                                    '--driver',stage/'tests/krylov_policy_driver','--output',root/'policy-timing'],stage)
        run(profile,'diagnostics',[sys.executable,'-B',stage/'tests/diagnostics_reproducibility.py',
                                   '--driver',stage/'tests/krylov_diagnostics_driver','--before-driver',before_driver,
                                   '--exe',stage/'source/QS3ED2','--output',root/'diagnostics','--performance'],stage)
        for label,cmd in [
            ('lapack',['make','check-lapack']),('standard',['make','check']),
            ('stability',['make','check-lanczos','LANCZOS_OUTPUT='+str(root/'stability')]),
            ('fields',['make','check-sparse-field','SPARSE_FIELD_OUTPUT='+str(root/'fields')]),
            ('stress',['make','check-krylov-stress','KRYLOV_OUTPUT='+str(root/'stress')]),
            ('precision',['make','check-precision-contract','PRECISION_OUTPUT='+str(root/'precision')]),
            ('independent',[sys.executable,'-B','tests/numerical_baseline.py','--exe',stage/'source/QS3ED2','--output',root/'independent']),
            ('extended',['make','check-long','top_srcdir='+str(stage),'top_builddir='+str(stage)]),
            ('historical',[sys.executable,'-B','tests/compare_historical.py','--build',stage,'--output',root/'historical.json'])]:
            if a.quick and label in ['extended','historical']:continue
            run(profile,label,cmd,stage)
        if a.quick:
            save();continue
        # Existing seeded CLI suite, instrumenting only an external exact source
        # snapshot. The ordinary helper's historical default revisions are intact.
        sys.path.insert(0,str(stage/'tests'))
        import krylov_benchmarks as kb
        kb.REVISIONS={'phase1c':BASE};probe=root/'probe';probe.mkdir()
        try:
            bins=kb.build(repo,probe,fc,flags+' '+a.ldflags,source_archives={'phase1c':after_archive})
            (probe/'provenance.json').write_text(json.dumps(dict(base=BASE,source_archive_sha256=sha(after_archive),
                                                               compiler=fc,flags=flags,ldflags=a.ldflags,
                                                               executable_sha256=sha(bins['phase1c'])),indent=2)+'\n')
            run(profile,'acceptance',[sys.executable,'-B',stage/'tests/krylov_acceptance_probe.py',
                                      '--exe',bins['phase1c'],'--output',root/'acceptance'],stage)
        except subprocess.CalledProcessError as error:
            manifest['commands'].append(dict(profile=profile,label='probe-build',exit_code=error.returncode,required=True))
            save()
        finally:
            # Avoid reusing profile-specific imports when a different compiler
            # profile starts. The helper itself contains no compiler state.
            sys.path.pop(0)
        save()
    failures=[r for r in manifest['commands'] if r['exit_code']!=0 and r['required']]
    print('Required command failures:',len(failures),flush=True)
    return int(bool(failures) or len(manifest['profiles'])!=len(a.profiles))


if __name__=='__main__':raise SystemExit(main())
