#!/usr/bin/env python3
"""Phase 1-E: characterize precision contracts without changing either solver.

Exit zero means the experiment and its diagnostic assertions completed, NOT that
all solves converged or met the proposed accuracy goals. report.json retains each
failure/warning. Requires only the existing optional NumPy test dependency.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import numpy as np
import krylov_stress as stress
from numerical_baseline import model
from sparse_field_validation import BONDS, FIXTURES


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical():
    fields = np.zeros((4, 3))
    for site, field in FIXTURES['single-site4']:
        fields[site] += field
    h, spin = model(4, BONDS, fields)
    e, u = np.linalg.eigh(h)
    return h, e, u, spin


def diagnostics(h, e, u, q, values):
    """Independent reference assessment; reference gaps are NOT online estimates."""
    residuals = np.linalg.norm(h @ q - q * values, axis=0)
    scale = float(max(abs(e)))
    width = float((e[-1] - e[0]) / 2)
    wanted = q.shape[1]
    # A numerical cluster definition for this experiment, not a production policy.
    cluster_eps = 64 * np.finfo(float).eps * max(width, np.finfo(float).tiny)
    members = np.flatnonzero(e <= e[wanted - 1] + cluster_eps)
    v = u[:, members]
    leakage = float(np.linalg.norm(q - v @ (v.conj().T @ q), 2))
    boundary_gap = float(e[wanted] - e[wanted - 1]) if wanted < len(e) else None
    complete = len(members) == wanted
    projector = float(np.linalg.norm(q @ q.conj().T - v @ v.conj().T, 2)) if complete else None
    # All references are known: distance from Ritz values to the unwanted spectrum.
    outside = np.setdiff1d(np.arange(len(e)), members)
    separation = float(np.min(abs(values[:, None] - e[outside]))) if len(outside) else None
    reference_defect = h @ u - u * e
    reference_error = float(np.linalg.norm(reference_defect, 2))
    # Generated analytic eigensystems also have formation/roundoff error. Include
    # their unwanted-space action defect; this is still a floating assessment,
    # not a validated interval enclosure of the exact real-arithmetic spectrum.
    reference_correction = float(np.linalg.norm(reference_defect[:, outside].conj().T @ q)) if len(outside) else 0.
    # Frobenius block-residual bound is safe here; do not assume spectral-norm
    # residual / gap universally bounds the spectral-norm Sylvester solution.
    bound = min(1., (float(np.linalg.norm(residuals)) + reference_correction) / separation) if separation and separation > 0 else None
    broad = u[:, e <= e[wanted - 1] + 1e-8 * max(width, np.finfo(float).tiny)]
    result = dict(residual_norms=residuals.tolist(), spectral_scale=scale,
                  centered_scale=width, boundary_gap=boundary_gap,
                  reference_separation=separation, reference_subspace_bound=bound,
                  reference_eigensystem_residual=reference_error,
                  reference_bound_correction=reference_correction,
                  orthogonality_error=float(np.linalg.norm(q.conj().T @ q - np.eye(wanted), 2)),
                  requested_space_leakage=leakage, projector_error=projector,
                  complete_cluster=complete,
                  broad_cluster_leakage=float(np.linalg.norm(q - broad @ (broad.conj().T @ q), 2)),
                  lowest_eigenvalue_errors=abs(values - e[:wanted]).tolist(),
                  absolute_residual_le_1e_minus10=bool(max(residuals) <= 1e-10),
                  relative_residual_le_1e_minus10=bool(max(residuals) <= 1e-10 * max(scale, np.finfo(float).tiny)),
                  centered_residual_le_1e_minus10=bool(max(residuals) <= 1e-10 * max(width, np.finfo(float).tiny)),
                  reference_subspace_le_1e_minus5=bool(complete and bound is not None and bound <= 1e-5))
    if bound is not None and complete:
        # Allow serialization/matmul roundoff, not a relaxed solver assertion.
        if leakage > bound + 2e-11:
            raise ValueError('independent subspace bound violated')
    return result


def observables(q, spin):
    psi = q[:, 0]
    lm = np.array([[np.vdot(psi, a @ psi).real for a in site] for site in spin])
    cf = np.array([[np.vdot(psi, spin[i][a] @ spin[j][b] @ psi).real
                    for b in range(3) for a in range(3)]
                   for i in range(4) for j in range(i + 1, 4)])
    return lm, cf


def analytic_controls():
    """Known states, explicitly not outputs of a changed/prospective solver."""
    e = np.array([-1., -1. + 1e-10, 0., 2.])
    h = np.diag(e); u = np.eye(4, dtype=complex)
    q = np.array([1., 1., 0., 0.], complex)[:, None] / np.sqrt(2.)
    value = np.real(np.diag(q.conj().T @ h @ q))
    mixed = diagnostics(h, e, u, q, value)
    mixed['energy_change_for_identical_candidates'] = 0.
    assert mixed['absolute_residual_le_1e_minus10'] and not mixed['reference_subspace_le_1e_minus5']
    assert mixed['projector_error'] > .7
    shifted = diagnostics(h + 1e6 * np.eye(4), e + 1e6, u, q, value + 1e6)
    # Exact excited state: every residual threshold accepts it; inclusion fails.
    excluded = diagnostics(h, e, u, u[:, 1:2], e[1:2])
    assert excluded['residual_norms'] == [0.] and excluded['projector_error'] == 1.
    # Degenerate ground cluster: full projectors agree under unitary rotations.
    deg_e = np.array([-1., -1., 0., 2.]); deg_h = np.diag(deg_e)
    rotation = np.array([[1., 1.], [-1., 1.]]) / np.sqrt(2.)
    full = diagnostics(deg_h, deg_e, u, u[:, :2] @ rotation, deg_e[:2])
    partial = diagnostics(deg_h, deg_e, u, q, deg_e[:1])
    assert full['projector_error'] < 1e-14 and partial['projector_error'] is None
    h, e, u, spin = physical()
    a = 5e-8
    q = (np.sqrt(1. - a*a) * u[:, 0] + a * u[:, 1])[:, None]
    value = np.real(np.diag(q.conj().T @ h @ q))
    state = diagnostics(h, e, u, q, value)
    lm, cf = observables(q, spin); exact_lm, exact_cf = observables(u[:, :1], spin)
    state.update(magnetization_error=float(np.max(abs(lm - exact_lm))),
                 correlation_error=float(np.max(abs(cf - exact_cf))),
                 legacy_residual_envelope_at_tol_1e14=1e-7 * max(abs(e)))
    assert max(state['residual_norms']) < state['legacy_residual_envelope_at_tol_1e14']
    assert not np.allclose(lm, exact_lm, atol=1e-9, rtol=1e-9)
    assert max(state['lowest_eigenvalue_errors']) < 1e-13
    return dict(cluster_mixture=mixed, shifted_cluster_mixture=shifted,
                excluded_ground=excluded, degenerate_complete=full,
                degenerate_partial=partial, constructed_physical_state=state)


def definitions():
    n = 32
    rng = np.random.default_rng(8223)
    u, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)))
    for gap in [1e-2, 1e-6, 1e-10, 0.]:
        e = np.r_[-1. + gap * np.arange(4), np.linspace(-.5, 2., n - 4)]
        h = (u * e) @ u.conj().T
        for alg in [1, 2]:
            for tol in [1e-12, 1e-20]:
                for seed in [19, 73]:
                    yield dict(family=f'gap-{gap}', n=n, h=h, e=e, u=u, alg=alg,
                               tol=tol, seed=seed, wanted=1 if alg == 1 else 4,
                               keep=8, capacity=16, budget=n if alg == 1 else 120)
    for family, n, h, e, u in stress.fixtures():
        if n != 32 or family != 'separated':
            continue
        for factor, shift in [(0., 0.), (0., 2.), (1e-9, 0.), (1e9, 0.), (1., 1e6)]:
            for alg in [1, 2]:
                yield dict(family=f'scale-{factor}-shift-{shift}', n=n,
                           h=factor*h+shift*np.eye(n), e=factor*e+shift, u=u,
                           alg=alg, tol=1e-12, seed=19, wanted=1 if alg == 1 else 3,
                           keep=6, capacity=12, budget=n if alg == 1 else 120)
    for item in stress.definitions():
        if item['family'] in ['excluded-ground', 'weak-ground-overlap']:
            yield item
    # Same independently constructed physical matrix as the GNU 13 failure.
    h, e, u, _ = physical()
    for alg in [1, 2]:
        for tol in [1e-14, 1e-18]:
            for seed in [19, 73, 211]:
                yield dict(family='physical-random', n=16, h=h, e=e, u=u, alg=alg,
                           tol=tol, seed=seed, wanted=1 if alg == 1 else 2,
                           keep=4, capacity=8, budget=16 if alg == 1 else 1000)
    # Explicit starts retain a small ground component while sampling the
    # orthogonal complement. This is a kernel experiment, not the old CLI RNG.
    for amplitude in [1e-2, 1e-4, 1e-6, 1e-8]:
        for seed in [19, 73, 211]:
            rng = np.random.default_rng(seed)
            coeff = rng.normal(size=15) + 1j*rng.normal(size=15)
            coeff /= np.linalg.norm(coeff)
            start = amplitude*u[:, 0]+np.sqrt(1.-amplitude**2)*(u[:, 1:]@coeff)
            yield dict(family=f'physical-weak-{amplitude}', n=16, h=h, e=e, u=u,
                       alg=1, tol=1e-14, seed=seed, wanted=1, keep=4, capacity=8,
                       budget=16, start=start)


def assess(driver, root, item, threads, index):
    work = root/f'{index:04d}-omp{threads}'; work.mkdir()
    wanted = item['wanted']; n = item['n']
    row = {k: v for k, v in item.items() if k not in ['h', 'e', 'u', 'start']}
    row.update(id=work.name, threads=threads)
    start = item.get('start', np.zeros(n, complex))
    header = (f"&problem n={n},alg={item['alg']},wanted={wanted},keep={item['keep']},"
              f"capacity={item['capacity']},maxiter={item['budget']},miniter=1,interval=1,"
              f"tol={item['tol']},seed={item['seed']},use_start={'.true.' if 'start' in item else '.false.'}/\n")
    deck = header + ''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in item['h'].flatten())
    deck += ''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in start)
    (work/'input.dat').write_text(deck)
    row['input_sha256'] = sha(work/'input.dat')
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (work/'input.dat').open('rb') as inp, (work/'stdout.log').open('wb') as log:
        proc = subprocess.run([str(driver)], cwd=work, env=env, stdin=inp,
                              stdout=log, stderr=subprocess.STDOUT, timeout=90)
    lines = (work/'status.dat').read_text().splitlines()
    status = int(lines[0])
    if status not in [-1, 0, 1, 2]:
        raise ValueError('unknown kernel status')
    if (status == 0) != (proc.returncode == 0):
        raise ValueError('driver exit and kernel status disagree')
    row.update(solver_status=status, reason=lines[1], iterations=int(lines[2]),
               exit_code=proc.returncode, legacy_success=status == 0,
               rng_state=(work/'rng.dat').read_text().splitlines(),
               h_applications=int((work/'metrics.dat').read_text().split()[0]))
    row['classification'] = {1: 'nonconvergence', 2: 'numerical_failure', -1: 'invalid_parameters'}.get(status, 'success')
    if status == 0:
        values = np.loadtxt(work/'eigenpairs.dat', ndmin=2)[:, 1]
        data = np.loadtxt(work/'vectors.dat', ndmin=2)
        q = (data[:, 2] + 1j*data[:, 3]).reshape(wanted, n).T
    elif status == 1 and item['alg'] == 2:
        values = np.loadtxt(work/'last-estimates.dat', ndmin=2)[:, 0]
        data = np.loadtxt(work/'last-vectors.dat', ndmin=2)
        q = (data[:, 0] + 1j*data[:, 1]).reshape(wanted, n).T
    else:
        row['candidate_assessment'] = None
        return row
    assessment = diagnostics(item['h'], item['e'], item['u'], q, values)
    if not np.isfinite(q).all() or not np.isfinite(values).all():
        raise ValueError('nonfinite candidate')
    if status == 0:
        row['classification'] = ('accepted_reference_targets' if
                                 max(assessment['lowest_eigenvalue_errors']) <= 1e-10 * max(assessment['spectral_scale'], 1.) and
                                 assessment['requested_space_leakage'] < 1e-5 else
                                 'accepted_accuracy_or_inclusion_warning')
    if item['family'].startswith('physical'):
        _, _, _, spin = physical()
        lm, cf = observables(q, spin); exact_lm, exact_cf = observables(item['u'][:, :1], spin)
        ok = np.allclose(lm, exact_lm, atol=1e-9, rtol=1e-9) and np.allclose(cf, exact_cf, atol=1e-9, rtol=1e-9)
        assessment.update(magnetization_error=float(np.max(abs(lm-exact_lm))),
                          correlation_error=float(np.max(abs(cf-exact_cf))),
                          unchanged_observable_assertions_met=bool(ok))
        sep = float(np.min(abs(values[0] - item['e'][1:])))
        angle_bound = min(1., assessment['residual_norms'][0]/sep) if sep > 0 else None
        # For a normalized pure state, |delta <A>| <= 2 ||A|| sin(theta).
        # ||S_i^a||=1/2 and ||S_i^a S_j^b||=1/4 for these distinct sites.
        assessment.update(reference_ground_separation=sep,
                          reference_magnetization_error_bound=angle_bound,
                          reference_correlation_error_bound=.5*angle_bound if angle_bound is not None else None)
    row['candidate_assessment'] = assessment
    # Compare criteria on the SAME endpoint. This does not simulate an earlier
    # stop, a new recurrence, or the cost of a prospective online policy.
    row['legacy_success'] = status == 0
    trace = re.findall(r'itr, ene, conv:\s+\d+\s+([\d.eE+-]+)\s+([\d.eE+-]+)', (work/'stdout.log').read_text())
    row['last_printed_energy_change'] = float(trace[-1][1]) if trace else None
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--driver', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    driver = args.driver.resolve(strict=True)
    root = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='ed2-precision-'))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    controls = analytic_controls()
    rows = []
    for index, item in enumerate(definitions()):
        for threads in [1, 2]:
            try:
                rows.append(assess(driver, root, item, threads, index))
            except (OSError, ValueError, AssertionError, subprocess.TimeoutExpired) as error:
                rows.append(dict(id=f'{index:04d}-omp{threads}', family=item['family'],
                                 classification='harness_or_runtime_error', error=str(error)))
    counts = dict(Counter(row['classification'] for row in rows))
    physical_failures = [r['id'] for r in rows if r.get('legacy_success') and
                         r['family'].startswith('physical') and
                         not r['candidate_assessment']['unchanged_observable_assertions_met']]
    report = dict(schema=1, interpretation='characterization; zero exit is not universal solver accuracy',
                  numpy=np.__version__, driver_sha256=sha(driver), harness_sha256=sha(Path(__file__)),
                  support_sha256={name: sha(Path(__file__).with_name(name)) for name in
                                  ['krylov_stress.py', 'numerical_baseline.py', 'sparse_field_validation.py']},
                  analytic_controls=controls, counts=counts,
                  successful_physical_accuracy_warnings=physical_failures, cases=rows)
    (root/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(total=len(rows), counts=counts, physical_accuracy_warnings=physical_failures)))
    return int(counts.get('harness_or_runtime_error', 0) > 0)


if __name__ == '__main__':
    raise SystemExit(main())
