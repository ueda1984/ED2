# 1000-site cubic O_h cache-capacity experiment (2026-11-18)

## Motivation and numerical provenance

PRESTO job `19383.jpmjpr1911`, pinned original
`7d36eacfa91b3c4a954a83b2c07c8d6e62c5c7da`, was independently
audited (single pass, operator attestation pending):

- 10x10x10 3D periodic J1-J2 XXZ spin 1/2, exactly 3 down spins
- Gamma / cubic O_h A1g, D = 4299
- E = 847.3078940655731; direct Ritz residual 2.2401986516701506e-9
- H applications = 170
- Maximum RSS = 6,053,960 KiB (5.773506 GiB)
- basis = 12.543675 s, operator_cache = 67.866437 s,
  lanczos = 0.124927 s, total = 80.557621 s.

The original cache builder `source/ham2vec.f90::make_wk_loc_and_ele`
uses `wk_dim=lv`, `MNTE=-1`, and allocates two dense rectangular arrays
`wk_loc(wk_dim1,wk_dim)` (default integer, 4 bytes) and
`wk_ele(wk_dim1,wk_dim)` (complex(8), 16 bytes) with
`wk_dim1=2*NO_one + 8*NO_two + 1`.

There are `NO_one=0`, `NO_two=9*1000=9000` bonds, and
`wk_dim=4299`. Thus
`(72001)*(4299)*(4+16)=6,190,645,980` bytes
(`5.765488 GiB`) are reserved and initialized for these
two arrays alone: **99.861% of observed RSS**. The majority of slots
are zero and never used in this fixed-magnetization model.

## Safe cap for the specified physical model

All terms conserve the number of down spins (Jx=Jy, D=0,
no transverse fields). Every down spin has at most six nearest and
twelve face-diagonal next-nearest neighbors. With exactly 3 down
spins there are **at most 3*(6+12)=54 directed down/up flip-flops**
per spin configuration (double-occupied edges are excluded).
The cache also stores one diagonal entry. The existing routine stores
the diagonal in slot 1, begins non-diagonal writes at 2, and allocates
`wk_dim1=MNTE+1`. Therefore **MNTE=64** leaves eleven extra
entries beyond the maximum 55 occupied cache slots, and is safe for
this exact benchmark input. It is **not** a safe universal bound
for arbitrary ED2 Hamiltonians, transverse fields, larger magnetization
sectors or altered lattice degree.

Expected cache array capacity with cap:
`(64+1)*4299*20=5,588,700` bytes = `5.3298 MiB`,
about **1/1107.7** the original reserved cache bytes.
This is an **analytical storage estimate**, not a measured total-RSS
or time claim. The code still loops over all 9000 bonds for each
of the 4299 representatives, so lower memory consumption need not
improve the 67.9 s operator construction time proportionally.

## Scope of code changes

This draft PR makes **no Fortran source changes**. In
`tests/cpc_cubic_pointgroup.py`, the generated input accepts
an optional `mnte` (default `-1`); the benchmark has
`--mnte` (default `-1`) and the guarded PBS launcher
reads `ED2_MNTE` only when explicitly supplied.

An independent audit re-reads `input.dat` and insists its
`MNTE` value equals the recorded trial and summary metadata;
the complete input bytes, source revision and executable are
still separately hashed. Command line validation enforces
`--mnte -1` or `--mnte >=54` for this specific three-magnon cube
model. The original `-1` behavior remains unchanged.

CI directly compares both `MNTE=-1` and `MNTE=64` using
the **same Fortran binary** and matching N=216 A1g/A2g
input model, and asserts agreement of reduced dimension, eigenenergy
(1e-8), certified direct residual and H-application count.
It logs both RSS and cache wall time as observations, without
a machine-dependent performance gate.

## PRESTO controlled experiment: N=1000, A1g, 1 thread

After CI and the new staged source SHA are recorded, and **only with
operator-approved exclusive-host use**, perform one pilot with
`ED2_MNTE=64` and all other physical parameters equal to job
`19383.jpmjpr1911`:

```sh
qsub -v "ED2_STAGE=$HOME/ED2_CUBIC_OH_mnte_stage,ED2_RESULTS_ROOT=$HOME/ED2_CUBIC_OH_results,ED2_PRESET=cubic-main,ED2_SIZES=10,ED2_THREADS=1,ED2_SECTORS=oh-a1g,ED2_WARMUPS=0,ED2_REPEATS=1,ED2_TIMEOUT=1800,ED2_MNTE=64,ED2_EXCLUSIVE_CONFIRMED=1" \
  tools/presto_cubic_exclusive.pbs
```

An independent auditor must pass and verify same energy (within
1e-8), residual (at most 1e-8), D=4299 and same H applications
before interpreting RSS/time. Verify the changed input's
`MNTE=64` is recorded, and separately preserve the original
19383 output. A one-off measurement is a **mechanistic pilot**,
not publication-grade scaling evidence.

## Next bottleneck if timing remains high

`make_wk_loc_and_ele` traverses all `D*NO_two =
4299*9000=38,691,000` physical edges per construction.
A future optional local-neighbor specialized builder could inspect
only occupied-spin neighbors (at most `4299*3*18 = 232,146`
candidate directed spin flips). The diagonal ZZ contribution
can be computed exactly as the polarized vacuum energy
minus per-down-site weighted degree/2 plus the ZZ couplings
between pairs of down spins. This **must** be independently
crosschecked at the Hamiltonian matrix level and with
fast/legacy symmetry and multiple small/large sectors before
integration. No such local-neighbor refactor is included here.

Do not claim a speedup until *measured* across repetitions,
with identical physical inputs, comparable core affinity and verified
exclusive-node allocation. Do not merge the experiment to main
without independent scientific and code review.
