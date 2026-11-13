
# Provisional example: response to symmetry-breaking fields

!!! warning "Historical analysis: missing evidence"
    The originally referenced figures `cutoff_test.png`, `fitting_functions.png`, and
    `finite_size_scaling.png` are not present in the repository.
    The numerical estimates and scaling claims below have **not** been
    independently reproduced in Phase 1 validation. This page is a
    historical methodological description, **not** a verified performance
    benchmark or a citable numerical result.

This historical note sketches a procedure for analyzing symmetry-breaking responses in quantum spin systems near saturation magnetization. It should not be treated as evidence of reproduced numerical results until the original data and scripts are restored.

The intended purpose is methodological: to show which calculations a user would need to perform with QS³-ED2, without asserting that the quoted numerical estimates have been independently validated in the current code revision.

In principle, such a workflow would use ED2 for

- calculations close to **saturation magnetization**
- controlled truncation of the Hilbert space using the **number of down spins**
- systematic analysis of **responses to small transverse fields**
- **finite‑size scaling analysis** of susceptibilities and order parameters

for lattice sizes that are difficult to treat using conventional exact‑diagonalization approaches.

The historical example concerns the **antiferromagnetic Heisenberg model on the square lattice** near the saturation field.

---

# Model

We consider the spin‑1/2 Heisenberg Hamiltonian

$$
\hat H =
J\sum_{\langle i,j\rangle}
\mathbf S_i \cdot \mathbf S_j
-
H \sum_i S_i^z
-
\delta h_u \sum_i S_i^x
-
\sum_i \delta h_s^{(i)} S_i^x .
$$

Here

- $J>0$ is the exchange interaction,
- $H$ is the longitudinal magnetic field,
- $\delta h_u$ is a **uniform transverse field**, and
- $\delta h_s^{(i)}$ is a **staggered transverse field**.

The transverse fields are introduced as **fictitious symmetry-breaking fields**
to probe the response of the system.

!!! note "ED2 field convention"
    This displayed model uses physical Zeeman minus signs. The ED2
    `FILE_hvec` implementation instead contributes `+ h · S`.
    Supply `h_z=-H` and the corresponding negative transverse
    field coefficients to represent the model above.

All calculations are performed **near the saturation magnetization**, where the number of down spins is small and QS³‑ED2 is particularly efficient.

---

# Hilbert‑Space Truncation Using the Number of Down Spins

QS³‑ED2 constructs the Hilbert space by specifying a cutoff

$$
N_{\downarrow}
$$

for the number of down spins.

Because the system is studied close to the saturation magnetization, the relevant states are contained in sectors with **small $N_{\downarrow}$**.

The first numerical task is therefore to determine **how large the cutoff $N_{\downarrow}$ must be** when small transverse perturbations are applied.

*The historical cutoff-convergence plot is not available in this repository.*

The historical analysis considered the transverse magnetization response to small **uniform** and **staggered** fields at

$$
H/J = 3.7 .
$$

The earlier notes treat this as a representative field within a finite-size magnetization interval; the precise plateau interpretation requires the missing reference data.

The finite-size selection of the longitudinal field, and any thermodynamic extrapolation of that value, have not been independently verified here.

The historical notes suggest testing an excitation cutoff extending to

$$
N_{\downarrow}^{\prime} + 1
$$

where $N_{\downarrow}^{\prime}$ is the number of down spins in the ground state when the transverse fields are zero.

This is physically reasonable because the transverse perturbation only mixes the ground state with nearby sectors.

Another important observation is that the response to the **staggered transverse field** is much larger than the response to the uniform transverse field.  
This behavior is physically natural because the staggered field couples directly to the spin pattern associated with the symmetry‑broken state.

---

# Extraction of the Order Parameter

The second step is to estimate the **spontaneous staggered magnetization** in the thermodynamic limit.

Because exact diagonalization can only be performed on finite systems, we analyze the dependence of the transverse magnetization on the small staggered field.

*The historical fit-comparison plot is not available in this repository.*

Two fitting functions are introduced:

$$
y_1 = a\tanh(g\,\delta h_s/J) + b(\delta h_s/J)
$$

$$
y_2 = y_1 + c(\delta h_s/J)^3 .
$$

These functions do not have a strict theoretical justification; they are introduced as **empirical fitting forms** to test the stability of the extrapolation.

Interpreting the fit parameter $a$ as a spontaneous order parameter
requires verifying the order of limits: take the thermodynamic limit
$L\to\infty$ **before** the vanishing staggered-field limit
$\delta h_s\to 0^+$. For a finite $g$, the stated hyperbolic-tangent
function vanishes as the field goes to zero; $a$ is therefore **not
automatically** the zero-field finite-size magnetization.

These empirical fitting functions are illustrative and cannot establish
the existence or value of spontaneous order without reliable size
scaling, fit-window tests and uncertainty estimates.

---

# Finite‑Size Scaling

The final step is to analyze the **finite‑size scaling** of the numerical results.

*The historical finite-size-scaling plot is not available in this repository.*

Two quantities are examined:

1. **Staggered susceptibility**. With sublattice pattern
   $\eta_i=\pm 1$ and $\delta h_s^{(i)}=\eta_i\delta h_s$,
   define $m_s^x=(1/N)\sum_i\eta_i\langle S_i^x\rangle$.
   Its zero-field linear response is

$
\chi_s=
\left.\frac{\partial m_s^x}{\partial \delta h_s}\right|_{\delta h_s=0}.
$

A finite difference $m_s^x/\delta h_s$ is only appropriate if the
zero-field baseline is zero and the applied field is sufficiently small.
The staggered sign factor must not be omitted.

2. **Spontaneous staggered magnetization** extracted from the fitting parameter $a$.

The unverified historical draft suggests an approximate scaling

$$
\chi_s \sim L^4 ,
$$

within the system sizes studied here. At present, we regard this as an empirical scaling observed in the present numerical data. A detailed comparison with previous theoretical and numerical studies remains to be carried out.

The historical draft quotes an extrapolated order parameter near

$$
m_x \approx 0.19
$$

at

$$
H/J = 3.79 \pm 0.03 .
$$

This value appears to be compatible with previously reported numerical estimates, although a detailed literature comparison is left for future work.

---

# Significance for QS³‑ED2

This historical draft outlines a possible numerical procedure for studying symmetry-breaking responses with ED2; it does **not** establish that the quoted calculations have been reproduced with the current software.

In particular, QS³‑ED2 enables

- efficient calculations **near saturation magnetization**
- controlled truncation using the **number of down spins**
- numerical evaluation of **response functions**
- **finite-size scaling procedures**, subject to independent validation

for system sizes that are difficult to treat using conventional exact‑diagonalization approaches.

To turn this outline into a verified benchmark, restore the original data, scripts and plots, or regenerate them from a pinned ED2 revision and report convergence and fit uncertainties.
