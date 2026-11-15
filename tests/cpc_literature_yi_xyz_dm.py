#!/usr/bin/env python3
"""Hamiltonian-level reproduction of Yi et al., PRB 100, 024423 (2019), Eq. (1).

This is a *finite-size, independent numerical cross-check*, not a
reproduction of published thermodynamic phase boundaries or WYSI results.
The paper uses Pauli sigma, periodic boundaries and a *minus* field term.
ED2 uses spin S=sigma/2 and a *plus* field term. A directed closing bond
(N,1) checks ED2's DM orientation and correlation-column conventions.

The main ED2 science demonstrator remains the alternating mixed-spin,
non-U(1) XYZ/Gamma/field/chirality response. This complementary test is
not a replacement by a simpler spin-conserving model.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

import numpy as np

# Written independently from the ED2 input couplings.
PAULI=(
    np.array([[0,1],[1,0]],complex),
    np.array([[0,-1j],[1j,0]],complex),
    np.diag([1.,-1.]).astype(complex),
)
J=1.0
GAMMA=0.37
DELTA=0.41
FIELD=0.58
D_VALUES=(0.,0.19)
ENERGY_TOL=3e-9
OBS_TOL=3e-9
FD_STEP=1e-4


def full_operators(n,local):
    identity=np.eye(2,dtype=complex)
    sites=[]
    for i in range(n):
        site=[]
        for op in local:
            m=np.array([[1.]],complex)
            for j in range(n):
                m=np.kron(m,op if i==j else identity)
            site.append(m)
        sites.append(site)
    return sites


def paper_hamiltonian(n,D):
    """Literal Eq. (1): J [(1+gamma)/2 sigx sigx + ...]
    + D (sigx sigy - sigy sigx) - h sigz, periodic oriented bonds.
    """
    ops=full_operators(n,PAULI)
    H=np.zeros((2**n,2**n),complex)
    K=np.zeros_like(H)
    pairs=[(i,(i+1)%n) for i in range(n)]
    for i in range(n):
        H-=FIELD*ops[i][2]
    for i,j in pairs:
        H+=J*((1+GAMMA)/2*(ops[i][0]@ops[j][0])
             +(1-GAMMA)/2*(ops[i][1]@ops[j][1])
             +DELTA*(ops[i][2]@ops[j][2]))
        K+=ops[i][0]@ops[j][1]-ops[i][1]@ops[j][0]
    H+=D*K
    return H,K,ops,pairs


def converted_ed2_data(n,D):
    """Exact coefficient conversion for S=sigma/2 and +hvec dot S."""
    Jx=2*J*(1+GAMMA)
    Jy=2*J*(1-GAMMA)
    Jz=4*J*DELTA
    bonds=[(i,(i+1)%n,[Jx,Jy,Jz,0.,0.,4*D,0.,0.,0.])
           for i in range(n)]
    fields=[[0.,0.,-2*FIELD] for _ in range(n)]
    return bonds,fields


def check_mapping(n,D,H):
    """Independently compare the paper's Pauli and the converted S forms."""
    spinops=full_operators(n,[x/2 for x in PAULI])
    bonds,fields=converted_ed2_data(n,D)
    mapped=np.zeros_like(H)
    for i in range(n):
        mapped+=fields[i][2]*spinops[i][2]
    for i,j,c in bonds:
        mapped+=(c[0]*spinops[i][0]@spinops[j][0]
                 +c[1]*spinops[i][1]@spinops[j][1]
                 +c[2]*spinops[i][2]@spinops[j][2]
                 +c[5]*(spinops[i][0]@spinops[j][1]
                        -spinops[i][1]@spinops[j][0]))
    error=float(np.max(np.abs(mapped-H)))
    if error>2e-12:
        raise AssertionError("Pauli-to-spin coefficient mapping error: "+str(error))
    return error


def reference(n,D):
    H,K,pauli,pairs=paper_hamiltonian(n,D)
    herm=float(np.max(np.abs(H-H.conj().T)))
    if herm>1e-12:
        raise AssertionError("Paper reference is non-Hermitian")
    mapping_error=check_mapping(n,D,H)
    vals,vecs=np.linalg.eigh(H)
    psi=vecs[:,0]
    spinops=[[x/2 for x in site] for site in pauli]
    def expect(A):
        z=np.vdot(psi,A@psi)
        if abs(z.imag)>1e-9:
            raise AssertionError("Complex expectation of Hermitian operator")
        return float(z.real)
    mag=np.array([[expect(spinops[i][a]) for a in range(3)]
                  for i in range(n)])
    corr=np.array([[expect(spinops[i][a]@spinops[j][b])
                    for b in range(3) for a in range(3)] for i,j in pairs])
    sz=sum((site[2] for site in spinops),np.zeros_like(H))
    comm=float(np.linalg.norm(H@sz-sz@H)/np.linalg.norm(H))
    if comm<1e-4:
        raise AssertionError("Anisotropic XYZ case unexpectedly U(1)-conserving")
    chirality=expect(K)
    Hplus=H+FD_STEP*K
    Hminus=H-FD_STEP*K
    hf_fd=float((np.linalg.eigvalsh(Hplus)[0]-
                 np.linalg.eigvalsh(Hminus)[0])/(2*FD_STEP))
    hf_error=abs(hf_fd-chirality)
    if hf_error>2e-5:
        raise AssertionError("Paper DM Hellmann-Feynman check failed")
    return dict(n=n,D=D,dim=2**n,bonds=converted_ed2_data(n,D)[0],
                fields=converted_ed2_data(n,D)[1],pairs=pairs,
                E=vals,mag=mag,corr=corr,gap=float(vals[1]-vals[0]),
                paper_K=chirality,HF_error=hf_error,
                nonU1_relative_commutator=comm,mapping_error=mapping_error)


def run_ed2(exe,root,ref,omp):
    n,D=ref["n"],ref["D"]
    tag=f"yi-prb2019-pbc-n{n}-D{D:.2f}-omp{omp}"
    work=root/tag
    work.mkdir(exist_ok=False)
    bonds,fields=ref["bonds"],ref["fields"]
    (work/"bonds.dat").write_text("".join(
        f"{i+1} {j+1} "+" ".join(f"{x:.17g}" for x in c)+"\n"
        for i,j,c in bonds))
    (work/"fields.dat").write_text("".join(
        f"{i+1} "+" ".join(f"{x:.17g}" for x in v)+"\n"
        for i,v in enumerate(fields)))
    (work/"spins.dat").write_text("0.5\n"*n)
    (work/"caps.dat").write_text("1\n"*n)
    # Keep the wrap-around pair in its original (N,1) order. The program
    # internally swaps it and must restore directed Cartesian output.
    (work/"pairs.dat").write_text("".join(
        f"{i+1} {j+1}\n" for i,j in ref["pairs"]))
    deck=f"""&input_parameters
NOS={n}, NODmin=0, NODmax={n}, ALG=3,
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={n}, NO_two={len(bonds)}, wk_dim=-1, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={n}, FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-12, minitr=5, maxitr=1000, itrint=10
/
&input_TRLan
NOE=2, NOK=4, NOM=16, maxitr=1000, lnc_ene_conv=1d-12,
i_vec_min=1, i_vec_max=2
/
"""
    (work/"input.dat").write_text(deck)
    env=dict(os.environ,OMP_NUM_THREADS=str(omp),OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as inp:
        proc=subprocess.run([str(exe)],stdin=inp,cwd=work,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=240)
    (work/"stdout.log").write_bytes(proc.stdout)
    if proc.returncode:
        raise AssertionError("ED2 dense solver failed: "+
                             proc.stdout.decode(errors="replace")[-2500:])
    energies=np.loadtxt(work/"eigenvalues.dat",ndmin=2)
    mags=np.loadtxt(work/"local_mag.dat",ndmin=2)
    corr=np.loadtxt(work/"two_body_cf_xyz.dat",ndmin=2)
    if energies.shape!=(ref["dim"],2) or mags.shape!=(n,4) or corr.shape!=(n,11):
        raise AssertionError("Unexpected ED2 dense output dimensions: "+
                             str((energies.shape,mags.shape,corr.shape)))
    np.testing.assert_array_equal(energies[:,0],np.arange(1,ref["dim"]+1))
    np.testing.assert_array_equal(mags[:,0],np.arange(1,n+1))
    np.testing.assert_array_equal(corr[:,:2],
                                  np.array([(i+1,j+1) for i,j in ref["pairs"]]))
    energy_err=float(np.max(np.abs(energies[:,1]-ref["E"])))
    mag_err=float(np.max(np.abs(mags[:,1:4]-ref["mag"])))
    corr_err=float(np.max(np.abs(corr[:,2:11]-ref["corr"])))
    # ED2 Cartesian file order: site_i,site_j,xx,yx,zx,xy,yy,zy,xz,yz,zz.
    # K in the paper uses Pauli; multiply ED2 xy-yx by four.
    ed2_K=4*float(np.sum(corr[:,5]-corr[:,3]))
    k_error=abs(ed2_K-ref["paper_K"])
    if max(energy_err,mag_err,corr_err,k_error)>max(ENERGY_TOL,OBS_TOL):
        raise AssertionError("Published-model independent mismatch: "+
                             str((energy_err,mag_err,corr_err,k_error)))
    return dict(case=tag,n_sites=n,dimension=ref["dim"],D=D,omp=omp,
                first_gap=ref["gap"],E0=float(ref["E"][0]),
                max_spectrum_error=energy_err,
                max_magnetization_error=mag_err,
                max_cartesian_correlation_error=corr_err,
                paper_chirality_error=k_error,
                HF_error=ref["HF_error"],
                nonU1_relative_commutator=ref["nonU1_relative_commutator"],
                pauli_spin_mapping_error=ref["mapping_error"])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="4,6,8")
    p.add_argument("--omp",type=int,choices=(1,2),default=1)
    a=p.parse_args()
    sizes=[int(x.strip()) for x in a.sizes.split(",")]
    if not sizes or any(x not in (4,6,8) for x in sizes) or len(set(sizes))!=len(sizes):
        raise ValueError("Allowed unique sizes: 4,6,8")
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(exist_ok=False,parents=True)
    rows=[]
    failures=[]
    for n in sizes:
        for D in D_VALUES:
            try:
                r=reference(n,D)
                rows.append(run_ed2(exe,root,r,a.omp))
                print("LITERATURE VALIDATED",rows[-1],flush=True)
            except Exception as exc:
                failures.append(dict(n=n,D=D,error=str(exc)))
                print("LITERATURE FAILURE",failures[-1],flush=True)
    report=dict(schema=1,source=os.environ.get("GITHUB_SHA","record git revision"),
                article="Yi et al., Physical Review B 100, 024423 (2019), Eq. (1)",
                doi="10.1103/PhysRevB.100.024423",
                scope="Hamiltonian and finite-size spectrum/observables only",
                parameters=dict(J=J,gamma=GAMMA,Delta=DELTA,h=FIELD,
                                D_values=D_VALUES,PBC=True,pauli_convention=True),
                successful=len(rows),expected=2*len(sizes),
                cases=rows,violations=failures)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("LITERATURE CHECK",len(rows),"/",2*len(sizes),"violations",len(failures),
          flush=True)
    return 0 if len(rows)==2*len(sizes) and not failures else 1


if __name__=="__main__":
    raise SystemExit(main())
