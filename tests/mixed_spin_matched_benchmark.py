#!/usr/bin/env python3
"""Matched-target Krylov characterization of the ED2 mixed-spin model.

Both algorithms are asked for *one* low eigenpair, use an identical
explicit complex start in ED2's own solver basis, and the same requested
absolute residual threshold. ALG2 may additionally generate seeded random
directions. The program does not silently promote legitimate nonconvergence
to an accepted physical result.

This is exploratory evidence: no changes to ED2 numerical kernels.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from mixed_spin_nonu1_demo import couplings
from mixed_spin_sparse_scale import reference_model

SEED=19
TOLS=(1e-6,1e-8,1e-10,1e-12)
FIELD=0.15
ENERGY_ATOL=3e-7
OBS_ATOL=8e-6


def diagnostic_fields(path):
    values={}
    for line in path.read_text().splitlines():
        if "=" in line:
            k,v=line.split("=",1)
            values[k.strip()]=v.strip()
    if "status" not in values or "h_applications" not in values:
        raise AssertionError("Missing required Krylov status and counters")
    return values


def initial_vector_text(dim):
    rng=np.random.default_rng(1984)
    re=rng.standard_normal(dim)
    im=rng.standard_normal(dim)
    return "".join(f"{a:.17g} {b:.17g}\n" for a,b in zip(re,im))


def setup_input(work,n,alg,model,initial_text):
    bonds,fields,spins,pairs=(model[x] for x in ("bonds","fields","spin","pairs"))
    (work/"bonds.dat").write_text("".join(
        f"{i+1} {j+1} "+" ".join(f"{x:.17g}" for x in c)+"\n"
        for i,j,c in bonds))
    (work/"fields.dat").write_text("".join(
        f"{i+1} "+" ".join(f"{x:.17g}" for x in vec)+"\n"
        for i,vec in enumerate(fields)))
    (work/"spins.dat").write_text("".join(f"{s:.1f}\n" for s in spins))
    (work/"caps.dat").write_text("".join(f"{round(2*s)}\n" for s in spins))
    (work/"pairs.dat").write_text("".join(f"{i+1} {j+1}\n" for i,j in pairs))
    (work/"initial-vector.txt").write_text(initial_text)
    nodmax=sum(round(2*s) for s in spins)
    deck=f"""&input_parameters
NOS={n}, NODmin=0, NODmax={nodmax}, ALG={alg},
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={n}, NO_two={len(bonds)}, wk_dim=-1, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={len(pairs)}, FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-12, minitr=5, maxitr=650, itrint=10
/
&input_TRLan
NOE=1, NOK=10, NOM=64, maxitr=260, lnc_ene_conv=1d-12,
i_vec_min=1, i_vec_max=1
/
"""
    (work/"input.dat").write_text(deck)


def run_one(exe,root,n,alg,eps,rep,model,initial,keep=10,capacity=64):
    tag=f"n{n}-alg{alg}-tol{eps:.0e}-repeat{rep}".replace("-","m",1)
    work=root/tag
    work.mkdir(exist_ok=False)
    setup_input(work,n,alg,model,initial)
    if alg==2:
        deck=(work/"input.dat").read_text()
        marker="NOE=1, NOK=10, NOM=64"
        if deck.count(marker)!=1:
            raise AssertionError("Ambiguous TRLan workspace template")
        deck=deck.replace(marker,f"NOE=1, NOK={keep}, NOM={capacity}")
        (work/"input.dat").write_text(deck)
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",
             MKL_NUM_THREADS="1",ED2_RNG_SEED=str(SEED),
             ED2_INITIAL_VECTOR=str(work/"initial-vector.txt"),
             ED2_DIAGNOSTICS="diagnostics.txt",ED2_RESIDUAL_ATOL=str(eps))
    timelog=work/"gnu-time.txt"
    cmd=[str(exe)]
    if Path("/usr/bin/time").exists():
        cmd=["/usr/bin/time","-f","MAX_RSS_KIB=%M","-o",str(timelog)]+cmd
    begin=time.perf_counter()
    with (work/"input.dat").open("rb") as stream:
        proc=subprocess.run(cmd,cwd=work,stdin=stream,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=260)
    elapsed=time.perf_counter()-begin
    (work/"stdout.log").write_bytes(proc.stdout)
    diagfile=work/"diagnostics.txt"
    if not diagfile.is_file():
        raise AssertionError("Missing diagnostic sidecar: "+proc.stdout.decode(errors="replace")[-1500:])
    diag=diagnostic_fields(diagfile)
    status=int(diag["status"])
    if (status==0)!=(proc.returncode==0):
        raise AssertionError("Process exit code does not match Krylov status")
    if diag.get("explicit_initial")!="T" or diag.get("seed_specified")!="T":
        raise AssertionError("Explicit start and seed not recorded")
    if int(diag["wanted"])!=1:
        raise AssertionError("Expected one requested eigenpair for both methods")
    if status not in (0,1,2):
        raise AssertionError("Unexpected solver status "+str(status))
    if alg==2 and (int(diag["keep"])!=keep or int(diag["capacity"])!=capacity):
        raise AssertionError("Krylov workspace not equal to requested values")
    results=dict(case=tag,n_sites=n,dimension=model["dimension"],
                 alg=alg,wanted=1,requested_keep=keep if alg==2 else None,
                 requested_capacity=capacity if alg==2 else None,
                 initial_sha256=hashlib.sha256(initial.encode()).hexdigest(),
                 seed=SEED,residual_atol=eps,repeat=rep,
                 status=status,reason=diag["reason"],exit_code=proc.returncode,
                 h_applications=int(diag["h_applications"]),
                 iterations=int(diag["iterations"]),cycles=int(diag["cycles"]),
                 wall_seconds=elapsed,peak_rss_kib=None,
                 lowest_state_inclusion=diag.get("lowest_state_inclusion"))
    if timelog.is_file():
        for line in timelog.read_text().splitlines():
            if line.startswith("MAX_RSS_KIB="):
                results["peak_rss_kib"]=int(line.split("=",1)[1])
    if status!=0:
        # The solver's failure is a *recorded numerical observation*.
        # It must not have generated a successful eigenpair/observable file.
        if any((work/name).exists() for name in
               ("eigenvalues.dat","local_mag.dat","two_body_cf_xyz.dat")):
            raise AssertionError("Nonconvergent solver wrote accepted outputs")
        return results

    eigen=np.loadtxt(work/"eigenvalues.dat",ndmin=2)
    if eigen.shape!=(1,2) or eigen[0,0]!=1:
        raise AssertionError("Eigenpair output does not match one-state goal")
    E=float(eigen[0,1])
    Eerror=abs(E-model["energies"][0])
    if Eerror>ENERGY_ATOL:
        raise AssertionError("Accepted wrong lowest state despite small residual "+
                             str((E,Eerror,diag["reason"])))
    mag=np.loadtxt(work/"local_mag.dat",ndmin=2)
    cf=np.loadtxt(work/"two_body_cf_xyz.dat",ndmin=2)
    if mag.shape!=(n,4) or cf.shape!=(len(model["pairs"]),11):
        raise AssertionError("Physical observables have wrong shape")
    np.testing.assert_array_equal(mag[:,0],np.arange(1,n+1))
    np.testing.assert_array_equal(cf[:,:2],
                                  np.array([(i+1,j+1) for i,j in model["pairs"]]))
    magerr=float(np.max(np.abs(mag[:,1:]-model["mag"])))
    cferr=float(np.max(np.abs(cf[:,2:]-model["corr"])))
    if max(magerr,cferr)>OBS_ATOL:
        raise AssertionError("Accepted inaccurate observables: "+str((magerr,cferr)))
    if diag.get("all_requested_residuals_met")!="T" or diag.get("requested_residual_met_1")!="T":
        raise AssertionError("Accepted without meeting residual contract")
    # The algorithm's Ritz/physical residual appears as the third
    # floating column in the one-index pair diagnostics.
    pair_lines=[line.split() for line in diagfile.read_text().splitlines()
                if len(line.split())==8 and line.split()[0]=="1"]
    if len(pair_lines)!=1:
        raise AssertionError("Exactly one diagnostic pair required")
    residual=float(pair_lines[0][5])
    if residual>eps*(1+2e-14):
        raise AssertionError("Accepted residual exceeds requested threshold")
    results.update(energy=E,energy_error=Eerror,residual=residual,
                   magnetization_max_error=magerr,correlation_max_error=cferr)
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--sizes",default="10,12")
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    sizes=[int(x) for x in args.sizes.split(",")]
    if sizes not in ([10],[12],[10,12]):
        raise ValueError("Supported --sizes combinations: 10, 12, 10,12")
    runs=[]
    references=[]
    violations=[]
    for n in sizes:
        print(f"REFERENCE n={n}...",flush=True)
        model=reference_model(n,FIELD)
        references.append(dict(n_sites=n,dimension=model["dimension"],
                               eigenvalues=model["energies"].tolist(),
                               reference_residuals=model["eigenpair_residuals"]))
        initial=initial_vector_text(model["dimension"])
        grid=TOLS if n==10 else (1e-8,1e-10)
        repeats=2 if n==10 else 1
        for eps in grid:
            for alg in (1,2):
                for rep in range(1,repeats+1):
                    try:
                        entry=run_one(exe,root,n,alg,eps,rep,model,initial)
                        runs.append(entry)
                        print(f"MEASURED n={n} alg={alg} atol={eps:.0e} rep={rep}: "
                              f"status={entry['status']} h_calls={entry['h_applications']} "
                              f"r={entry.get('residual')} wall={entry['wall_seconds']:.3f}",
                              flush=True)
                    except Exception as error:
                        violations.append(dict(n_sites=n,alg=alg,atol=eps,rep=rep,error=str(error)))
                        print("ASSERTION FAILURE:",violations[-1],flush=True)
    report=dict(schema=1,kind="one-pair matched-start tolerance study",
                source_sha=os.environ.get("GITHUB_SHA","record in CI"),
                models=references,runs=runs,contract_violations=violations,
                run_count=len(runs),accepted=sum(x["status"]==0 for x in runs),
                nonconverged=sum(x["status"]!=0 for x in runs))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("MATCHED BENCHMARK:",report["run_count"],"runs;",
          report["accepted"],"accepted;",report["nonconverged"],"nonconverged;",
          len(violations),"contract violations",flush=True)
    expected=sum(4*2*2 if n==10 else 2*2 for n in sizes)
    if len(runs)!=expected or violations:
        return 1
    return 0


if __name__=="__main__":
    raise SystemExit(main())
