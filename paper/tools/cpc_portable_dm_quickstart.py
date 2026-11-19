#!/usr/bin/env python3
"""Portable ED2 CPiP example: finite-size DM/chirality response, stdlib only.

Runs the *connected*, two-site-periodic, alternating spin-1/2/spin-1
non-U(1) XYZ/Gamma/staggered-DM Hamiltonian at D=(-h,0,+h).
Every ED2 job uses the complete Hilbert space (NODmin=0 and
NODmax=sum_i(2S_i)); there is no unphysical fixed-Sz restriction.

For D=0 the ED2 ground energy, oriented Cartesian correlation
chirality K = sum_j (-1)^j (xy-yx), its central energy derivative,
and finite-size susceptibility are compared to archival independently
constructed sparse/dense reference results.

NO NumPy, SciPy, Matplotlib or Python wheel installation is needed
for this *end-user* reproduction (unlike the independent validation
suite). The external executable needs its own Fortran BLAS/LAPACK runtime.

This reproduces a previously verified illustrative physics test, not
a new phase, a thermodynamic extrapolation, or a universal solver
certification. It deliberately rejects swapped xy/yx file columns.
"""
import argparse
import hashlib
import json
import math
import os
import re
from pathlib import Path
import subprocess
import sys
import time

DEFAULT_STEP=2e-4
REFERENCE={
    4: dict(E0=-2.2309372506,K=0.001150247666,chi=0.6343112017685585),
    6: dict(E0=-3.58386349587,K=0.001090053482,chi=0.5687329253277751),
}
COLS=("xx","yx","zx","xy","yy","zy","xz","yz","zz")
SEED=19
RESIDUAL_REQUEST=1e-8


def to_float(raw):
    x=float(raw.replace("D","E").replace("d","e"))
    if not math.isfinite(x):raise ValueError("Nonfinite numeric output")
    return x


def model(size,D):
    if size not in REFERENCE:raise ValueError("Validated standalone sizes are 4 and 6")
    spins=[.5,1.]*(size//2)
    bonds=[]
    for j in range(size-1):
        sign=(-1)**j
        coeffs=[.77+.05*sign,.63+.03*sign,1.09-.04*sign,
                0.,0.,D*sign,0.,0.,.065+.015*sign]
        bonds.append((j+1,j+2,coeffs))
    fields=[]
    for i in range(size):
        sign=(-1)**i
        fields.append((i+1,[.15*(.4*sign+.1),0.,.17+.035*sign]))
    return spins,bonds,fields


def write_input(work,size,D):
    spins,bonds,fields=model(size,D)
    work.mkdir(parents=True,exist_ok=False)
    (work/"spins.dat").write_text(
        "".join(f"{v:.1f}\n" for v in spins),encoding="utf-8")
    (work/"caps.dat").write_text(
        "".join(f"{int(round(2*v))}\n" for v in spins),encoding="utf-8")
    (work/"bonds.dat").write_text(
        "".join(f"{i} {j} "+" ".join(f"{x:.17g}" for x in vals)+"\n"
                for i,j,vals in bonds),encoding="utf-8")
    (work/"fields.dat").write_text(
        "".join(f"{i} "+" ".join(f"{x:.17g}" for x in vals)+"\n"
                for i,vals in fields),encoding="utf-8")
    (work/"pairs.dat").write_text(
        "".join(f"{j} {j+1}\n" for j in range(1,size)),encoding="utf-8")
    caps=sum(int(round(2*v)) for v in spins)
    deck=f"""&input_parameters
NOS={size}, NODmin=0, NODmax={caps}, ALG=2,
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={size}, NO_two={size-1}, wk_dim=-1, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={size-1}, FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-12, minitr=5, maxitr=1200, itrint=10
/
&input_TRLan
NOE=2, NOK=4, NOM=24, maxitr=320, lnc_ene_conv=1d-12,
i_vec_min=1, i_vec_max=2
/
"""
    (work/"input.dat").write_text(deck,encoding="utf-8")


def read_numeric(path,columns,rows):
    if not path.is_file():
        raise AssertionError("ED2 output missing: "+str(path))
    values=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():continue
        fields=line.split()
        if len(fields)!=columns:
            raise AssertionError(f"Unexpected {path.name} column count: {len(fields)}")
        values.append([to_float(x) for x in fields])
    if len(values)!=rows:
        raise AssertionError(f"Unexpected {path.name} row count: {len(values)}")
    return values


def diagnostics(path):
    if not path.is_file():
        raise AssertionError("Missing opt-in solver diagnostic sidecar")
    scalar={}
    pairs=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k,v=line.split("=",1)
            scalar[k.strip()]=v.strip()
        else:
            fields=line.split()
            if len(fields)==8 and fields[0].isdigit():
                pairs.append(dict(index=int(fields[0]),available=fields[1],
                                  certified=fields[2],ritz=to_float(fields[3]),
                                  residual=to_float(fields[5])))
    return scalar,pairs



def legacy_main_solver_diagnostics(stdout,dimension,energies):
    """Original QS3ED2 TRLan has no direct-residual diagnostic sidecar.

    Accept only a converged restart indicator AND two independently
    computed native <Phi|H|Phi> expectations consistent with Ritz output.
    This is NOT a direct ||H psi-E psi|| residual certification.
    The separate reference E0/K/chi checks remain mandatory.
    """
    if "Eigen solver: my thick restart Lanczos" not in stdout:
        raise AssertionError("Expected original-main thick-restart solver")
    dim=re.findall(r"(?m)^\s*THS\(k\)\s*=\s*(\d+)\s*$",stdout)
    if len(dim)!=1 or int(dim[0])!=dimension:
        raise AssertionError("Incorrect native projected Hilbert dimension")
    restart=re.findall(
        r"(?m)^tritr, ene\(NOK\), max_conv, flag:\s*"
        r"(\d+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+(\d+)\s*$",stdout)
    if not restart:
        raise AssertionError("Missing original-main TRLan restart diagnostic")
    iteration,_,metric,flag=restart[-1]
    conv=to_float(metric)
    if int(flag)!=2 or not (0.0<=conv<=1e-10):
        raise AssertionError("TRLan did not report converged restart: "+
                             str(restart[-1]))
    marker="Check accuracy <Phi|H|Phi>"
    if stdout.count(marker)!=1:
        raise AssertionError("Native eigenvector expectation check missing")
    tail=stdout.split(marker,1)[1].split("total lanczos step:",1)[0]
    rows=re.findall(r"(?m)^\s*([12])\s+\(\s*"
                    r"([-+0-9.EeDd]+)\s*,\s*([-+0-9.EeDd]+)\s*\)\s*$",tail)
    if len(rows)!=2 or [int(row[0]) for row in rows]!=[1,2]:
        raise AssertionError("Unexpected native H expectation table")
    checks=[complex(to_float(re),to_float(im)) for _,re,im in rows]
    discrepancy=max(abs(z.real-e) for z,e in zip(checks,energies))
    maximum_imag=max(abs(z.imag) for z in checks)
    if discrepancy>1e-8 or maximum_imag>1e-8:
        raise AssertionError("Native Ritz/energy expectation consistency failed: "+
                             str((discrepancy,maximum_imag)))
    steps=re.findall(r"(?m)^\s*total lanczos step:\s*(\d+)\s*$",stdout)
    if len(steps)!=1 or int(steps[0])<1:
        raise AssertionError("Missing original-main Lanczos iteration count")
    return dict(
        kind="original_main_no_direct_residual_sidecar",
        direct_residual_certified=False,
        max_direct_residual=None,
        native_restart_iteration=int(iteration),
        native_restart_metric=conv,
        max_ritz_expectation_difference=discrepancy,
        max_imaginary_ritz_expectation=maximum_imag,
        lanczos_steps=int(steps[0]),
        h_applications=None,
        seed_verified=False,
    )

def run(exe,root,size,D,threads,timeout):
    label=("zero" if D==0 else f"{D:+.10g}".replace("+","p").replace("-","m"))
    work=root/f"D-{label}"
    write_input(work,size,D)
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             ED2_RNG_SEED=str(SEED),
             ED2_RESIDUAL_ATOL=str(RESIDUAL_REQUEST),
             ED2_DIAGNOSTICS="diagnostics.txt")
    begin=time.perf_counter()
    with (work/"input.dat").open("rb") as source:
        proc=subprocess.run([str(exe)],cwd=work,stdin=source,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            check=False,timeout=timeout)
    elapsed=time.perf_counter()-begin
    (work/"stdout.log").write_bytes(proc.stdout)
    if proc.returncode!=0:
        raise AssertionError("Native ED2 exited with error code {}: {}".format(
            proc.returncode, proc.stdout.decode("utf-8",errors="replace")[-2000:]))
    dimension=6**(size//2)
    energy_rows=read_numeric(work/"eigenvalues.dat",2,2)
    if [int(v[0]) for v in energy_rows]!=[1,2]:
        raise AssertionError("Unexpected lowest eigenpair indices")
    energies=[v[1] for v in energy_rows]
    if (work/"diagnostics.txt").is_file():
        # Enhanced experimental branch: direct residuals are certified by
        # a source-produced diagnostic file, not deduced from energies.
        diag,pairs=diagnostics(work/"diagnostics.txt")
        status=int(diag["status"])
        if status!=0:
            raise AssertionError("ED2 failed convergence/exit contract: "+
                                 str((status,proc.returncode,diag.get("reason"))))
        if int(diag["dimension"])!=dimension or int(diag["wanted"])!=2:
            raise AssertionError("Incorrect full Hilbert dimension or requested pairs")
        if diag.get("all_requested_residuals_met")!="T" or \
           diag.get("residual_policy_enabled")!="T" or diag.get("seed_specified")!="T":
            raise AssertionError("Solver residual/seed policy not honored")
        if len(pairs)!=2 or any(p["index"]!=i+1 or p["available"]!="T" or
                                p["certified"]!="T" or
                                p["residual"]>RESIDUAL_REQUEST*(1+2e-12)
                                for i,p in enumerate(pairs)):
            raise AssertionError("Invalid low-eigenpair diagnostics")
        evidence=dict(kind="enhanced_direct_residual_sidecar",
                      direct_residual_certified=True,
                      max_direct_residual=max(p["residual"] for p in pairs),
                      h_applications=int(diag["h_applications"]),
                      seed_verified=True,
                      lanczos_steps=None)
    else:
        # The actual current generic main has no ED2_DIAGNOSTICS implementation.
        # Its restart metric and <Phi|H|Phi> agreement are useful checks,
        # but are *not* a direct eigenpair norm residual.
        evidence=legacy_main_solver_diagnostics(
            proc.stdout.decode("utf-8",errors="replace"),dimension,energies)
    status=0
    if energies[1]-energies[0]<.05:
        raise AssertionError("Unexpectedly small physical lowest gap")
    mag=read_numeric(work/"local_mag.dat",4,size)
    if [int(r[0]) for r in mag]!=list(range(1,size+1)):
        raise AssertionError("Incorrect site magnetization indices")
    correlations=read_numeric(work/"two_body_cf_xyz.dat",11,size-1)
    total=0.
    bond_values=[]
    for j,row in enumerate(correlations):
        if (int(row[0]),int(row[1]))!=(j+1,j+2):
            raise AssertionError("Wrong requested nearest-neighbor correlation pair")
        # ED2 files explicitly use xx,yx,zx,xy,yy,zy,xz,yz,zz:
        xy=row[5]
        yx=row[3]
        component=xy-yx
        bond_values.append(component)
        total+=((-1)**j)*component
    return dict(D=D,energy0=energies[0],first_gap=energies[1]-energies[0],
                kappa_sum=total,bond_kappa=bond_values,
                diagnostic=evidence,
                wall_seconds=elapsed,case=str(work.relative_to(root)),
                status=status,dimension=dimension)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path,
                   help="Existing, compiled ED2 QS3ED2 executable")
    p.add_argument("--output",required=True,type=Path,
                   help="New directory; keep complete input and outputs")
    p.add_argument("--size",type=int,choices=(4,6),default=4)
    p.add_argument("--step",type=float,default=DEFAULT_STEP)
    p.add_argument("--omp",type=int,choices=(1,2),default=1)
    p.add_argument("--timeout",type=int,default=180)
    args=p.parse_args()
    if not (math.isfinite(args.step) and abs(args.step-DEFAULT_STEP)<1e-16):
        p.error("The frozen independent-reference acceptance is valid only at step=0.0002")
    exe=args.exe.resolve(strict=True)
    if not exe.is_file():p.error("ED2 executable is not a file")
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    rows=[];issues=[]
    for D in (-args.step,0.,+args.step):
        try:
            result=run(exe,root,args.size,D,args.omp,args.timeout)
            rows.append(result)
            print("PORTABLE ED2 DM",result["case"],"E0",result["energy0"],
                  "K",result["kappa_sum"],"evidence",
                  result["diagnostic"]["kind"],flush=True)
        except Exception as exc:
            problem=dict(D=D,error_type=type(exc).__name__,message=str(exc))
            issues.append(problem)
            print("PORTABLE ED2 FAILURE",problem,flush=True)
    metrics={}
    if len(rows)==3:
        left,middle,right=rows
        h=args.step
        hf=(right["energy0"]-left["energy0"])/(2*h)
        chi=-(right["energy0"]-2*middle["energy0"]+left["energy0"])/(h*h*(args.size-1))
        ref=REFERENCE[args.size]
        metrics=dict(
            K_from_ED2_correlations=middle["kappa_sum"],
            central_dE_dD=hf,
            HF_abs_error=abs(hf-middle["kappa_sum"]),
            susceptibility_per_bond=chi,
            reference_E0=ref["E0"],
            reference_K=ref["K"],
            reference_chi=ref["chi"],
            E0_abs_error=abs(middle["energy0"]-ref["E0"]),
            K_abs_error=abs(middle["kappa_sum"]-ref["K"]),
            chi_abs_error=abs(chi-ref["chi"]),
        )
        if (metrics["HF_abs_error"]>3e-5 or
            metrics["E0_abs_error"]>1e-6 or
            metrics["K_abs_error"]>1e-5 or
            metrics["chi_abs_error"]>3e-3 or
            chi<0):
            issues.append(dict(error_type="PhysicsReferenceMismatch",
                               message="HF / independently archived E0, K or chi acceptance violated",
                               metrics=metrics))
    complete=len(rows)==3 and not issues
    def sha256(path):
        digest=hashlib.sha256()
        with path.open("rb") as src:
            for block in iter(lambda:src.read(1<<20),b""):
                digest.update(block)
        return digest.hexdigest()
    binary_hash=sha256(exe)
    report=dict(schema=2,source_sha=os.environ.get("GITHUB_SHA"),
                bin_path=str(exe),bin_sha256=binary_hash,
                model="connected two-site periodic non-U1 mixed-spin XYZ/Gamma/DM OBC",
                caveat="finite-size software demonstration, not a thermodynamic phase claim",
                size=args.size,omp=args.omp,finite_difference_step=args.step,
                spin_dimensions="alternating 2 and 3",independent_reference_source=
                "paper/data/cpc-stationary-dm-reference.csv (independent source; size-stationary rows)",
                rows=rows,metrics=metrics,failures=issues,
                pass_complete=complete)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("PORTABLE DM QUICKSTART",len(rows),"/3 jobs","pass",complete,
          "HF error",metrics.get("HF_abs_error"),
          "chi",metrics.get("susceptibility_per_bond"),flush=True)
    return 0 if complete else 1


if __name__=="__main__":
    raise SystemExit(main())
