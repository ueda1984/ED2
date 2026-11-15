#!/usr/bin/env python3
"""Second-order quasi-degenerate reference for the interacting eight-site ED2 stress model.

For the deliberately weak anisotropic XY bond
  epsilon*(S_edge^x S_core^x + eta*S_edge^y S_core^y)
we derive an effective 2x2 Hamiltonian within
span{|up> x |g_core>, |down> x |g_core>}.

Second-order self-energies are evaluated using TWO constrained sparse
linear solves in the orthogonal complement of the independently computed
seven-site core ground state. The augmented inverse is numerically
nonsingular even though H_core-E_core is singular.

The two-level model is not exact. It explicitly omits virtual-transition
terms of order O(epsilon^3/gap_core^2) and small delta-dependent corrections
to the denominators. It provides a distinct *algorithmic* reference to
the full 1296-dimensional sparse eigsh problem, NOT an independently
specified Hamiltonian or a certified high-precision exact answer.

This module is used only for the ultraweak, non-U(1), OBC mixed-spin
adversarial test and is not a general ED2 feature.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, splu
from mixed_spin_sparse_scale import construct_sparse_hamiltonian

ETA=0.6


def effective_second_order(model):
    """Compute low-doublet splitting and edge Sz through second order.

    Args:
        model: eight-site nonconserving sparse benchmark reference model,
               containing spins, bonds, fields, epsilon and delta.

    Returns:
        JSON-compatible scalar diagnostics. All input arrays are kept local.
    """
    eps=float(model["epsilon"])
    delta=float(model["delta"])
    if not (model["n_sites"]==8 and 0<eps<=1e-5 and abs(delta)<=1e-5):
        raise ValueError("Second-order reduction is supported only for the ultraweak eight-site probe")
    spins=model["spins"]
    bonds=model["bonds"]
    fields=model["fields"]
    if spins!=[0.5,1.0]*4 or len(bonds)!=7:
        raise AssertionError("Model does not have the required mixed-spin product structure")
    first=bonds[0]
    if first[:2]!=(0,1):
        raise AssertionError("Unexpected first bond")
    coeffs=np.asarray(first[2],float)
    if (abs(coeffs[0]-eps)>1e-14 or
        abs(coeffs[1]-ETA*eps)>1e-14 or
        np.max(np.abs(coeffs[2:]))>1e-14):
        raise AssertionError("This reduction assumes only epsilon*XX + eta*epsilon*YY on the edge")
    if not np.allclose(fields[0],[0.,0.,delta],rtol=0,atol=1e-15):
        raise AssertionError("Unexpected isolated-edge field")
    core_bonds=[(i-1,j-1,co) for i,j,co in bonds[1:]]
    if len(core_bonds)!=6 or any(i+1!=j for i,j,_ in core_bonds):
        raise AssertionError("Unexpected active-core open-chain graph")
    core,ops=construct_sparse_hamiltonian(spins[1:],core_bonds,fields[1:])
    d=core.shape[0]
    if d!=648:
        raise AssertionError("Unexpected active Hilbert dimension")
    eye=sparse.identity(d,dtype=complex,format="csc")
    rng=np.random.default_rng(191122)
    guess=rng.standard_normal(d)+1j*rng.standard_normal(d)
    guess/=np.linalg.norm(guess)
    energies,states=eigsh(core,k=2,which="SA",tol=1e-13,
                          ncv=min(90,d-1),maxiter=40000,v0=guess)
    order=np.argsort(energies)
    e0=float(energies[order[0]])
    e1=float(energies[order[1]])
    phi=states[:,order[0]]
    phi=phi/np.linalg.norm(phi)
    core_gap=e1-e0
    core_res=float(np.linalg.norm(core@phi-e0*phi))
    if not .2<core_gap<2. or core_res>5e-11:
        raise AssertionError("Active-core gap/residual is unsuitable for perturbative reduction")
    x=ops[0][0]
    y=ops[0][1]
    mx=float(np.vdot(phi,x@phi).real)
    my=float(np.vdot(phi,y@phi).real)
    if abs(np.vdot(phi,x@phi).imag)>1e-12 or abs(np.vdot(phi,y@phi).imag)>1e-12:
        raise AssertionError("Nonreal expectation of a Hermitian core operator")
    # Upper->lower transition from the sigma_x/y edge bond.
    # The upper state's outgoing Q excitation is (X+i*eta*Y)|g>.
    # The lower state's outgoing Q excitation is (X-i*eta*Y)|g>.
    up_source=(x@phi)+1j*ETA*(y@phi)
    dn_source=(x@phi)-1j*ETA*(y@phi)
    wup=up_source-phi*np.vdot(phi,up_source)
    wdn=dn_source-phi*np.vdot(phi,dn_source)

    # [Hc-e0  |g>][q] = [w], [<g| 0][lambda]=[0].
    # This solves q = Q (Hc-e0)^(-1) Q w without forming
    # a pseudoinverse or subtracting nearly degenerate full energies.
    col=sparse.csc_matrix(phi[:,None])
    row=sparse.csc_matrix(phi.conj()[None,:])
    aug=sparse.bmat([[core.tocsc()-e0*eye,col],
                     [row,sparse.csc_matrix((1,1),dtype=complex)]],format="csc")
    solver=splu(aug)

    def complement_solve(w):
        rhs=np.empty(d+1,dtype=complex)
        rhs[:d]=w
        rhs[d]=0
        sol=solver.solve(rhs)
        q=sol[:d]
        orth=abs(np.vdot(phi,q))
        projected_res=np.linalg.norm((core-e0*eye)@q-w)
        if orth>1e-9 or projected_res>2e-9:
            raise AssertionError("Constrained Q resolvent solve failed "+str((orth,projected_res)))
        expectation=float(np.vdot(w,q).real)
        if abs(np.vdot(w,q).imag)>1e-9 or expectation<=-1e-10:
            raise AssertionError("Nonpositive Q resolvent expectation")
        return q,expectation,float(projected_res),float(orth)

    qup,tup,rup,oup=complement_solve(wup)
    qdn,tdn,rdn,odn=complement_solve(wdn)
    # H_eff = [[delta/2 + Sigma_up, eps/2*(mx-i*eta*my)],
    #          [eps/2*(mx+i*eta*my), -delta/2+Sigma_down]]
    # where Sigma_up=-eps^2/4*<w_up|R|w_up>, and analogously down.
    shift_up=-(eps*eps/4.)*tup
    shift_down=-(eps*eps/4.)*tdn
    offdiag=(eps/2.)*(mx-1j*ETA*my)
    heff=np.array([[delta/2+shift_up,offdiag],
                   [offdiag.conjugate(),-delta/2+shift_down]],dtype=complex)
    if np.max(np.abs(heff-heff.conj().T))>5e-15:
        raise AssertionError("Effective matrix is non-Hermitian")
    e_eff,vec_eff=np.linalg.eigh(heff)
    gap_second=float(e_eff[1]-e_eff[0])
    gap_first=float(np.hypot(delta,eps*np.hypot(mx,ETA*my)))
    z_first=float(-delta/(2*gap_first))
    z_second_twolevel=float((abs(vec_eff[0,0])**2-
                              abs(vec_eff[1,0])**2)/2.)
    # First Q-space correction to the *state*: an up component
    # acquires a down/core-excited amplitude, and vice versa.
    up_amp=vec_eff[0,0]*phi-(eps/2.)*vec_eff[1,0]*qdn
    dn_amp=vec_eff[1,0]*phi-(eps/2.)*vec_eff[0,0]*qup
    probability_up=float(np.vdot(up_amp,up_amp).real)
    probability_down=float(np.vdot(dn_amp,dn_amp).real)
    z_second_dressed=(probability_up-probability_down)/(2*(probability_up+probability_down))
    if not 0<gap_second<.001 or not all(abs(z)<=.500000001 for z in
                       (z_first,z_second_twolevel,z_second_dressed)):
        raise AssertionError("Unphysical effective model values")
    if abs(z_second_twolevel-z_second_dressed)>2e-7:
        raise AssertionError("Dressed edge magnetization exceeds perturbative control")
    return dict(
        epsilon=eps,delta=delta,
        core_energy=e0,core_excitation_gap=core_gap,
        core_ground_direct_residual=core_res,
        core_mx=mx,core_my=my,
        sigma_up=shift_up,sigma_down=shift_down,
        sigma_difference=shift_up-shift_down,
        first_order_gap=gap_first,
        second_order_gap=gap_second,
        first_order_edge_z=z_first,
        second_order_edge_z=z_second_twolevel,
        second_order_dressed_edge_z=float(z_second_dressed),
        upper_Q_residual=rup,lower_Q_residual=rdn,
        upper_Q_orthogonality=oup,lower_Q_orthogonality=odn,
        upper_Q_norm=float(np.linalg.norm(qup)),
        lower_Q_norm=float(np.linalg.norm(qdn)),
        method="sparse constrained Q resolvent and 2x2 degenerate perturbation",
        validity="ultraweak eight-site modified-edge Hamiltonian only",
    )
