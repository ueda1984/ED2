# Input format

QS³-ED2 reads its input from **standard input (STDIN)** and expects a **Fortran NAMELIST** file (typically named `input.dat`).  
This page documents the **Fortran NAMELIST fields in ED2 source**.
Development-branch solver controls are described in
[Solver accuracy and reproducibility](solver-accuracy.md).

> **Reproducibility note**  
> Always record (i) the ED2 git tag/commit hash, (ii) `input.dat`, and (iii) all referenced `list_*.dat` files.

---

## Quick start

> **Important**
> The directories under `examples/` (e.g. `chain`, `cubic`, `kagome`, ...) are generated automatically
> after running `make check-long` on a fresh clone.  
> If they are missing, please run:
>
> ```bash
> make check-long
> ```
>
> before trying the commands below.

From an example directory:

```bash
cd examples/chain
../../source/QS3ED2 < input.dat > output.dat
```

The helpers under `script_*/run.sh` generate model inputs and run ED2.
For a fresh checkout, build ED2 first and execute the appropriate helper
to generate `examples/<model>/`; do not assume every generated example
directory already exists.

---

## Examples directory structure

The repository provides reference examples under `examples/`:

- `chain/`
- `cubic/`
- `cubic_sp_HB/`
- `honeycomb/`
- `j1j2_chain/`
- `kagome/`
- `mixed_spin_chain/`
- `square/`
- `triangular/`

After its model-generation helper has run, the example directory normally contains:

- `input.dat` (main NAMELIST input)
- model-definition lists such as `list_xyz_dm_gamma.dat`, `list_hvec.dat`, `list_spin.dat`, ...
- optional symmetry definition lists such as `list_p1.dat` ... `list_p6.dat`

We recommend starting from `examples/chain/` to understand the basic workflow.

---

## File format: Fortran NAMELIST

`input.dat` is a Fortran NAMELIST file consisting of multiple blocks:

- `&input_parameters ... &end` : system size, basis restriction, symmetry sector, I/O, etc.
- `&input_lancz ... &end`      : Lanczos parameters (used when `ALG=1`)
- `&input_TRLan ... &end`      : Thick-Restart Lanczos parameters (used when `ALG=2`)

### NAMELIST conventions (practical)

- Each block begins with `&<block_name>` and ends with `&end`.
- Entries are `KEY = VALUE` separated by commas (`,`). Newlines are allowed, but commas are recommended.
- Strings must be quoted, e.g., `"list_hvec.dat"` or `'list_hvec.dat'`.
- Double precision constants are written as `1.0d-14` (Fortran `d` exponent).
- Only keys declared in the source NAMELISTs are supported. An unknown name can cause a runtime input error; use the supplied examples and source declarations.

---

## Example: `input.dat` (from `examples/chain`)

```text
&input_parameters
  NOS = 100,
  NODmax = 3,
  NODmin = 0,

  L1 = 100,  L2 = 1,  L3 = 1,  L4 = 1,  L5 = 1,  L6 = 1,
  M1 = 0,    M2 = 0,  M3 = 0,  M4 = 0,  M5 = 0,  M6 = 0,

  NO_one = 100,
  NO_two = 100,

  wk_dim = -1,
  MNTE = -1,

  ALG = 1,
  wr_wf = 0,
  re_wf = 0,

  FILE_xyz_dm_gamma = "list_xyz_dm_gamma.dat",
  FILE_hvec         = "list_hvec.dat",

  FILE_wf = "./",
  OUTDIR  = "./",

  FILE_S1 = "list_p1.dat",
  FILE_S2 = "list_p2.dat",
  FILE_S3 = "list_p3.dat",
  FILE_S4 = "list_p4.dat",
  FILE_S5 = "list_p5.dat",
  FILE_S6 = "list_p6.dat",

  FILE_NODmax = "list_NODmax.dat",
  FILE_spin   = "list_spin.dat",

  cal_lm = 1,
  cal_cf = 1,

  NO_two_cf   = 9,
  FILE_two_cf = "./list_ij_cf.dat",
&end

&input_lancz
  lnc_ene_conv = 1.0d-14,
  minitr       = 5,
  maxitr       = 10000,
  itrint       = 5,
&end

&input_TRLan
  NOE          = 10,
  NOK          = 20,
  NOM          = 40,
  maxitr       = 10000,
  lnc_ene_conv = 1.0d-14,
  i_vec_min    = 1,
  i_vec_max    = 1,
&end
```

---

## `&input_parameters` reference

This block defines: (i) the target Hilbert subspace, (ii) optional symmetry sectors, (iii) Hamiltonian term lists, (iv) memory/performance knobs, and (v) output paths.

### System size and basis restriction

| Key | Type | Meaning |
|---|---:|---|
| `NOS` | INTEGER | Number of sites (spins), typically denoted by **N**. |
| `NODmin`, `NODmax` | INTEGER | Range of the number of spin-lowering operators, controlling the targeted subspace. See below. |

#### Basis convention (how the subspace is specified)

QS³-ED2 treats the fully polarized state $|v\rangle$ as a vacuum and spans the subspace by applying $N_{\downarrow}$ lowering operators:

$$
|a\rangle = \prod_{o=1}^{N_{\downarrow}} \hat{S}^{-}_{r_o}\,|v\rangle,\qquad
N_{\downarrow}\in[\mathrm{NODmin},\mathrm{NODmax}].
$$

- If you want the **full Hilbert space**, choose a wide range (model-dependent).  
- If your model has a **U(1) conservation law** (e.g., XXZ-like total $S^z$ conservation), using `NODmin=NODmax` selects a fixed sector (fixed total $S^z$).
- With excitation-nonconserving terms (e.g., transverse fields), NOD restriction defines a **projected Hamiltonian**, not generally an invariant sector of the unrestricted system.

Local constraints and local spin values are further specified by `FILE_NODmax` and `FILE_spin`.

---

### Symmetry sector specification (up to 6 commuting operations)

QS³-ED2 can incorporate up to **six commuting symmetry operations** (typically translations and point-group operations).  
You specify:

1) their **orders** `L1..L6` such that $\hat{T}_m^{L_m}=\hat{1}$  
2) their **sector labels** `M1..M6` such that

$$
\hat{T}_m|\Psi\rangle = e^{-i 2\pi M_m/L_m}\,|\Psi\rangle,\qquad M_m\in\{0,\dots,L_m-1\}.
$$

3) their **site permutations** via `FILE_S1..FILE_S6`

| Key | Type | Meaning |
|---|---:|---|
| `L1..L6` | INTEGER | Orders $L_m$. If `Lm=1`, the operator is treated as trivial and the corresponding `FILE_Sm` can be omitted. |
| `M1..M6` | INTEGER | Sector labels $M_m$ for each symmetry. |
| `FILE_S1..FILE_S6` | CHARACTER | Files encoding the site permutation induced by $\hat{T}_m$ (1-based indexing). |

**Example (1D translation by one site)**  
The permutation vector is $(2,3,\dots,N,1)$, stored as an `NOS`-line integer file.

> Tip  
> If you set a symmetry but the chosen sector is incompatible (empty), ED2 should report an error or produce no states in that sector.

---

### Hamiltonian term lists

!!! important "Field sign and bond orientation"
    A field entry contributes **+ h_x S^x + h_y S^y + h_z S^z** at
    the specified site. To reproduce a physical convention - B · S,
    supply h = -B.

    Each row of `FILE_hvec`: `site hx hy hz` (1-based site index).
    Each row of `FILE_xyz_dm_gamma`:
    `site_i site_j Jx Jy Jz Dx Dy Dz Gx Gy Gz`.
    The parser flips the three DM coefficients when reordering a
    reversed input bond. The field record number and its site label
    are different quantities.


| Key | Type | Meaning |
|---|---:|---|
| `NO_one` | INTEGER | Number of one-body terms (fields) read from `FILE_hvec`. |
| `NO_two` | INTEGER | Number of two-body terms (XYZ + DM + Γ, etc.) read from `FILE_xyz_dm_gamma`. |
| `FILE_hvec` | CHARACTER | Path to one-body term list. ED2 reads the first `NO_one` lines. |
| `FILE_xyz_dm_gamma` | CHARACTER | Path to two-body term list. ED2 reads the first `NO_two` lines. |

> Column formats depend on the model/example.  
> The safest reference is the corresponding `examples/<name>/list_*.dat` shipped in this repo.

---

### Memory/performance knobs for sparse Hamiltonian construction

| Key | Type | Meaning |
|---|---:|---|
| `wk_dim` | INTEGER | Controls how many sparse matrix elements are cached. `-1` enables automatic setting (recommended to start). |
| `MNTE` | INTEGER | Controls the internal limit on transitions per representative state. Too small values can trigger an error. `-1` enables automatic setting. |

**Practical recommendation:** start with `wk_dim = -1` and `MNTE = -1`.  
During runtime, ED2 prints an **Optimal MNTE**, which can be used to tune memory usage for large runs.

---

### Diagonalization algorithm and wavefunction I/O

| Key | Type | Meaning |
|---|---:|---|
| `ALG` | INTEGER | Algorithm selector: `1` Lanczos (ground state), `2` Thick-Restart Lanczos, `3` full diagonalization (LAPACK). |
| `wr_wf` | INTEGER | If `1`, write eigenvectors to `FILE_wf`. |
| `re_wf` | INTEGER | If `0`, diagonalize. If `1`, read eigenvectors from `FILE_wf` and compute observables only. |
| `FILE_wf` | CHARACTER | Directory for eigenvector I/O. |
| `OUTDIR` | CHARACTER | Output directory for observables and diagnostics. |

The Lanczos and Thick-Restart Lanczos algorithms used in QS³-ED2 follow standard formulations:

- C. Lanczos, *An iteration method for the solution of the eigenvalue problem of linear differential and integral operators*, J. Res. Natl. Bur. Stand. **45**, 255 (1950).
- K. Wu and H. Simon, *Thick-restart Lanczos method for large symmetric eigenvalue problems*, SIAM J. Matrix Anal. Appl. **22**, 602–616 (2000).

---

### Observable outputs (local magnetization, correlations)

| Key | Type | Meaning |
|---|---:|---|
| `cal_lm` | INTEGER | If `1`, compute local magnetization and write `OUTDIR/local_mag.dat`. |
| `cal_cf` | INTEGER | If `1`, compute two-point correlators and write `OUTDIR/two_body_cf_xyz.dat` and `OUTDIR/two_body_cf_z+-.dat`. |
| `NO_two_cf` | INTEGER | Number of site pairs $(r,r')$ for correlation functions. |
| `FILE_two_cf` | CHARACTER | File containing `NO_two_cf` pairs $(r,r')$. |

Typical output formats:

- `OUTDIR/local_mag.dat` stores $(r, \langle S^x_r\rangle, \langle S^y_r\rangle, \langle S^z_r\rangle)$.
- `OUTDIR/two_body_cf_xyz.dat` stores $(r, r', \langle S^\alpha_r S^\beta_{r'}\rangle)$ for $\alpha,\beta \in \{x,y,z\}$.
- `OUTDIR/two_body_cf_z+-.dat` stores correlators in the $(z,+,-)$ basis (see examples for exact column ordering).

---

## `FILE_NODmax` and `FILE_spin`

These two files define local constraints and local spin values.

### `FILE_NODmax` (local maximum number of lowering operators)

- Contains **`NOS` lines** of integers.
- Line `r` specifies the maximum number of lowering operators allowed on site `r`.
- Must satisfy:
  - not exceeding the global `NODmax`
  - not exceeding $2S_r$ on each site

Examples:
- Spin-1/2 chain: all lines are `1`
- Mixed-spin chain (e.g., alternating spins): pattern like `2,1,2,1,...`

### `FILE_spin` (local spin values)

- Contains **`NOS` lines** of real numbers.
- Line `r` gives the local spin $S_r$ (e.g., `0.5`, `1.0`, ...)

This design lets you target near-saturation sectors efficiently even for large-spin systems.

---

## `&input_lancz` (Lanczos)

This block is used when `ALG = 1`.

Common keys (as used in examples):

| Key | Type | Meaning |
|---|---:|---|
| `lnc_ene_conv` | REAL*8 | Energy convergence tolerance. |
| `minitr` | INTEGER | Minimum number of iterations. |
| `maxitr` | INTEGER | Maximum number of iterations. |
| `itrint` | INTEGER | Printing interval for iteration logs. |

---

## `&input_TRLan` (Thick-Restart Lanczos)

This block is used when `ALG = 2`.

Common keys (as used in examples):

| Key | Type | Meaning |
|---|---:|---|
| `NOE` | INTEGER | Number of eigenpairs to compute (requested). |
| `NOK`, `NOM` | INTEGER | TRLan workspace/restart parameters. |
| `maxitr` | INTEGER | Maximum number of iterations. |
| `lnc_ene_conv` | REAL*8 | Energy convergence tolerance. |
| `i_vec_min`, `i_vec_max` | INTEGER | Index range of vectors to output/use (example-dependent). |

---

## Best practices

1. Start from `examples/<model>/input.dat` and modify incrementally.
2. Keep referenced file paths (`FILE_*`) either:
   - relative to the directory where ED2 is run, or
   - absolute paths (for batch/HPC runs).
3. Record environment details (compiler, BLAS/LAPACK, thread settings) for reproducible comparisons.
4. See [Solver accuracy and reproducibility](solver-accuracy.md) before interpreting a solver's success as an observable-accuracy guarantee.

---

## Troubleshooting checklist

- **ED2 cannot find input/list files**  
  Confirm `FILE_*` paths are relative to the runtime working directory (or make them absolute).

- **Symmetry sector seems empty / ED2 stops early**  
  Verify `L*`, `M*` and the permutations in `FILE_S*` are consistent with your lattice.

- **Memory-related stop / slow runs**  
  Start with `wk_dim=-1`, `MNTE=-1`. Use runtime-reported `Optimal MNTE` to tune for large cases.
