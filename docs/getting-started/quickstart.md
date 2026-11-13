# Quickstart

This page demonstrates the **shortest build-and-smoke-test workflow** to verify that QS³-ED2
builds and runs correctly on your system.

The objective is deliberately minimal:

- build the code,
- execute a reference calculation,
- confirm successful output.

Conceptual background, design philosophy, and citation information are documented on the main documentation home page.

---

## 1. Clone and build

```bash
git clone https://github.com/QS-Cube/ED2.git
cd ED2
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j
```

This produces the executable:

```
source/QS3ED2
```

The binary is intentionally **not installed into PATH**.
All example scripts resolve it explicitly to avoid environment-dependent failures on HPC systems.

---

## 2. Run the smoke test

```bash
make check
```

This executes a minimal reference calculation and verifies basic functionality.

Expected result:

```
PASS: test_ex1.sh
```

A pass means this smoke-test path completed; it does not independently
validate every Hamiltonian, solver or physical observable.

---

## 3. Run a full example manually

Move to the chain example:

```bash
cd script_chain
./run.sh
```

This generates the corresponding chain input and performs an ED calculation, printing:

- solver configuration,
- convergence information,
- eigenvalues,
- local magnetization,
- selected correlation functions.

Output files such as:

```
eigenvalues.dat
local_mag.dat
two_body_cf_xyz.dat
```

are generated under `examples/chain/` by the helper script.

Important properties of the example scripts:

- no recompilation occurs,
- the executable is resolved explicitly,
- PATH is never used,
- helper programs are built at `make` time only.

These choices avoid path-related failures. Record the source revision, inputs, compiler, BLAS implementation, thread counts and initial-state settings when numerical reproducibility is required.

---

## 4. Optional: extended regression suite

```bash
make check-long
```

Runs all reference examples sequentially.

Recommended after:

- moving to a new machine or cluster,
- changing compiler,
- changing BLAS / MKL backend.

---

## 5. Thread control (recommended)

QS³-ED2 uses OpenMP in the main solver and threaded BLAS/LAPACK libraries.
To avoid nested parallelism during validation runs:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Users may later tune these values for performance benchmarking.
