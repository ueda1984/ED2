# ED2 (QS³-ED2): Exact Diagonalization for Quantum Spin Systems

**ED2 (QS³-ED2)** is a research-grade Fortran code for **Exact
Diagonalization of quantum spin lattice models**.

QS³ stands for **Quantum Spin System Solver**, reflecting the design
philosophy of ED2:

> **transparent physics, reproducible numerics, and controllable Hilbert
> spaces.**

ED2 is developed for:

-   **Reproducible computational physics**
-   **Explicit symmetry-resolved calculations**
-   **Controlled Hilbert-space truncation**
-   **High-performance shared-memory execution via OpenMP**

This documentation describes the ED2 source revision used to build the site.
Features in a development fork should not be assumed to exist in the
upstream QS-Cube/ED2 release until explicitly integrated.

Its primary objective is to guarantee that any third party can:

1.  Build ED2 from source,
2.  Run automated regression tests,
3.  Perform reference calculations,
4.  Understand all inputs and outputs,
5.  Reproduce published results,
6.  Cite the software unambiguously.

------------------------------------------------------------------------

## Overview

Exact Diagonalization remains one of the most reliable and transparent
approaches for strongly correlated quantum many-body systems.

ED2 focuses on **quantum spin lattice models** and supports both:

-   full Hilbert spaces, and\
-   **restricted (truncated) subspaces**,

allowing efficient calculations in dilute or symmetry-constrained
sectors.

Typical applications include:

-   Ground-state calculations\
-   Low-energy excitation spectra\
-   Local magnetization\
-   Two-body correlation functions

ED2 emphasizes **explicit construction of Hilbert spaces and
symmetries**, prioritizing clarity and reproducibility over black-box
automation.

------------------------------------------------------------------------

## Physical model and scientific scope

ED2 inherits the core design philosophy of the quantum spin solver QS3,
originally developed for studies of spin-1/2 XXZ systems near saturation
magnetization, and generalizes it to fully anisotropic quantum spin
Hamiltonians of the form

$$
\hat H =
\sum_{r>r'} \sum_{\alpha\in\{x,y,z\}}
\Big[
J^{\alpha}_{rr'} \hat S^{\alpha}_r \hat S^{\alpha}_{r'}
+ D^{\alpha}_{rr'} ( \hat{\mathbf S}_r \times \hat{\mathbf S}_{r'} )_{\alpha}
+ \Gamma^{\alpha}_{rr'} ( \hat S^{\beta}_r \hat S^{\gamma}_{r'} + \hat S^{\gamma}_r \hat S^{\beta}_{r'} )
\Big]
+ \sum_r \mathbf h_r \cdot \hat{\mathbf S}_r .
$$

with mutually distinct indices $\alpha\neq\beta\neq\gamma$. Here
$J^{\alpha}_{rr'}$, $D^{\alpha}_{rr'}$, and $\Gamma^{\alpha}_{rr'}$
denote exchange, Dzyaloshinskii--Moriya, and symmetric anisotropic
interactions between spins on sites $r,r'$, while $\mathbf h_r$
represents local Zeeman fields. Local spin magnitudes
$S_r = 1/2, 1, 3/2, \ldots$ may be specified independently at each site.

!!! important "Magnetic-field sign convention"
    ED2's field input enters the implemented Hamiltonian as **+ h · S**,
    not - h · S. To model a Zeeman term - B · S, supply h = -B
    in the field records (in your Hamiltonian's energy units).
    See [Input format](user-guide/input-format.md).

Typical scientific applications include:

-   finite-size scaling analyses of quantum phase transitions\
-   level spectroscopy of low-energy excitations\
-   identification of Berezinskii--Kosterlitz--Thouless transitions via
    twisted boundary conditions\
-   Anderson tower analyses of symmetry breaking

These techniques are widely used in numerical studies of low-dimensional
quantum magnetism.

------------------------------------------------------------------------

## Key features

### Explicit symmetry decomposition

Users explicitly specify symmetry generators (translations, rotations,
inversion, etc.) in the input file. ED2 constructs symmetry-resolved
Hilbert spaces, enabling block-diagonalization of the Hamiltonian and
substantial reductions in computational cost.

------------------------------------------------------------------------

### General quantum spin Hamiltonians

-   Heisenberg and XYZ-type interactions\
-   Dzyaloshinskii--Moriya and symmetric anisotropies\
-   External magnetic fields

------------------------------------------------------------------------

### Controlled Hilbert-space truncation

Calculations may be restricted to sectors with a fixed number of
spin-down excitations or equivalent constraints, dramatically reducing
memory and runtime.

------------------------------------------------------------------------

### Numerical solvers

-   Full diagonalization (small systems)
-   Lanczos and Thick-Restart Lanczos iterative solvers (larger systems)

------------------------------------------------------------------------

### Parallelization

Shared-memory parallelism via **OpenMP**.\
ED2 follows a single-node execution model (no MPI dependency).

------------------------------------------------------------------------

### Linear algebra backends

BLAS/LAPACK via:

-   OpenBLAS
-   Intel MKL
-   Vendor implementations

------------------------------------------------------------------------

## Typical workflow

``` bash
git clone https://github.com/QS-Cube/ED2.git
cd ED2
autoreconf -vfi
./configure FC=gfortran   # or FC=ifort / ifx
make -j
make check
make check-long
```

The commands

``` bash
make check
make check-long
```

provide automated regression tests and reference examples.

After successful compilation:

``` bash
cd examples/...
../source/QS3ED2 < input.dat > output.dat
```

Output files contain energies and physical observables.

------------------------------------------------------------------------

## Important threading note

ED2 uses OpenMP internally. When combined with threaded BLAS libraries
(e.g. MKL), nested parallelism should be avoided:

``` bash
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
```

Failure to do so may significantly degrade performance.

------------------------------------------------------------------------

## Reproducibility

For scientific publications, users are strongly encouraged to record:

-   ED2 version or Git commit hash
-   Compiler and version
-   BLAS/LAPACK backend
-   OpenMP settings
-   Input files

This information is essential for long-term reproducibility.

------------------------------------------------------------------------

## Citation

A journal manuscript describing ED2 is in preparation.

Until a DOI becomes available, please cite:

    ED2 (QS³-ED2), QS-Cube/ED2, GitHub repository,
    version <tag> or commit <hash>.

A `CITATION.cff` file is provided.

------------------------------------------------------------------------

## License

ED2 is released under the **MIT License**.

------------------------------------------------------------------------

## Authors and maintainers

-   Hiroshi Ueda\
-   Daisuke Yamamoto

------------------------------------------------------------------------

## Support and contributions

Bug reports and feature requests are welcome via GitHub Issues.
Contributions can be made through Pull Requests.

For scientific use, always accompany results with sufficient metadata to
guarantee reproducibility.
