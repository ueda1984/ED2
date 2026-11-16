#!/usr/bin/env python3
"""Independent 1000-site Gamma two-magnon Heisenberg reference.

Uses translation-orbit *relative displacement* states, NOT the
QS3/ED2 combinatorial state-list or Hamiltonian-vector code.
The cubic NN spin-1/2 Hamiltonian has Jxy=Jz=+1, 3000 periodic bonds.
The two-particle translation orbits are r equivalent to -r;
seven nonzero order-two displacements yield orbit size 500, while
the others have orbit size 1000. This gives 503 Gamma states.

Source-pinned original QS3 and ED2 both previously returned
E0 approximately 738.00677427221. We independently reconstruct
the 503x503 real symmetric matrix and verify it with SciPy eigsh.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

L=(10,10,10)
N=math.prod(L)
NEAREST=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
FROZEN_QS3_ED2_GROUND=738.00677427221

def add(a,b):
    return tuple((x+y)%m for x,y,m in zip(a,b,L))
def neg(a):
    return tuple((-x)%m for x,m in zip(a,L))
def canonical(a):
    return min(a,neg(a))
def orbsize(r):
    return N//2 if neg(r)==r else N

def construct():
    reps=sorted({canonical(r) for r in itertools.product(*(range(m) for m in L))
                 if r!=(0,0,0)})
    assert len(reps)==503
    assert sum(orbsize(r) for r in reps)==math.comb(N,2)
    assert sum(orbsize(r)==N//2 for r in reps)==7
    index={r:i for i,r in enumerate(reps)}
    H=sparse.lil_matrix((len(reps),len(reps)),dtype=np.float64)
    M=3*N
    for j,r in enumerate(reps):
        # H_ZZ = (M/4) - 2*z/2 + (# adjacent down-down bonds).
        # z=6; a nearest-neighbor pair restores one Ising exchange bond.
        H[j,j]=M/4.-6.+sum(r==step for step in NEAREST)
        for step in NEAREST:
            # Either flipped spin can move by one lattice bond.
            for r2 in (add(r,step),add(r,neg(step))):
                if r2==(0,0,0):
                    continue
                target=canonical(r2)
                H[index[target],j]+=0.5*math.sqrt(orbsize(r)/orbsize(target))
    H=H.tocsr()
    anti=H-H.transpose()
    violation=np.max(np.abs(anti.data),initial=0)
    if violation>1e-13:
        raise AssertionError(f"Independent matrix not Hermitian: {violation}")
    return H

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--report",type=Path,required=True)
    args=p.parse_args()
    H=construct()
    eigvals,eigvecs=eigsh(H,k=2,which="SA",tol=1e-13,maxiter=20000)
    order=np.argsort(eigvals)
    E0=float(eigvals[order[0]])
    E1=float(eigvals[order[1]])
    residual=float(np.linalg.norm(H@eigvecs[:,order[0]]-E0*eigvecs[:,order[0]]))
    err=abs(E0-FROZEN_QS3_ED2_GROUND)
    if residual>1e-9 or err>2e-8:
        raise AssertionError(f"Independent two-magnon crosscheck failed: residual={residual}, energy error={err}")
    out=dict(source_method="two-spin relative-displacement translation orbits",
             sites=N,shape=L,down_spins=2,k=(0,0,0),
             unsymmetrized_dimension=math.comb(N,2),
             reduced_dimension=H.shape[0],
             Hamiltonian_nnz=int(H.nnz),Jxy=1,Jz=1,
             first_energy=E0,second_energy=E1,gap=E1-E0,
             Ritz_residual=residual,
             source_pinned_old_new_energy=FROZEN_QS3_ED2_GROUND,
             independent_difference=err,
             note="Independent relative-separation exact Gamma Hamiltonian; "
                  "finite-size check not a speedup or bulk limit")
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print("QS3/ED2 INDEPENDENT 1000-SITE TWO-MAGNON",out,flush=True)

if __name__=="__main__":
    main()
