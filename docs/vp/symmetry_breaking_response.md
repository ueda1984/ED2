
# Provisional example: response to symmetry-breaking fields

!!! warning "Historical analysis: missing evidence"
    The originally referenced figures `cutoff_test.png`, `fitting_functions.png`, and
    `finite_size_scaling.png` are not present in the repository.
    The numerical estimates and scaling claims below have **not** been
    independently reproduced in Phase 1 validation. This page is a
    historical methodological description, **not** a verified performance
    benchmark or a citable numerical result.

This performance example demonstrates how **QS³‑ED2** can reproduce numerical analyses of symmetry‑breaking responses in quantum spin systems near saturation magnetization.

The purpose of this page is **not to present new physics**, but to demonstrate that the numerical procedure used in previous studies can be **practically executed using QS³‑ED2** for relatively large system sizes.

In particular, this benchmark illustrates that QS³‑ED2 enables

- calculations close to **saturation magnetization**
- controlled truncation of the Hilbert space using the **number of down spins**
- systematic analysis of **responses to small transverse fields**
- **finite‑size scaling analysis** of susceptibilities and order parameters

for lattice sizes that are difficult to treat using conventional exact‑diagonalization approaches.

The calculations below are performed for the **antiferromagnetic Heisenberg model on the square lattice** near the saturation field.

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

This value corresponds to the **center of the magnetization plateau** intersecting this value of $H/J$.

Using the plateau center stabilizes the numerical calculation when transverse fields are introduced.  
Although the optimal longitudinal field slightly shifts for finite systems, the extrapolated value **converges to $H/J = 3.7$ in the thermodynamic limit**, so the procedure remains consistent.

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

A key point is that in both cases the parameter

$$
a
$$

corresponds to the **staggered magnetization in the thermodynamic limit**

$$
m_x = a \qquad (L \rightarrow \infty).
$$

Thus, by fitting the finite‑field response we obtain an estimate of the spontaneous order parameter.

---

# Finite‑Size Scaling

The final step is to analyze the **finite‑size scaling** of the numerical results.

*The historical finite-size-scaling plot is not available in this repository.*

Two quantities are examined:

1. **Staggered susceptibility**

$$
\chi_s =
\frac{1}{N\,\delta h_s}
\sum_i \langle S_i^x \rangle
$$

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
