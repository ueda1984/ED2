# Direct complex CSR staging batch tuning: explicit memory/setup-time tradeoff

**Status:** small opt-in refinement of the already merged general complex
CSR two-pass builder, **not** a new Hamiltonian algorithm or a default
behavior change. All cubic-specific, Lanczos and CPC experimental
branches remain separate.

## Use

```bash
ED2_CSR_COMPLEX_CACHE=1 \
ED2_CSR_RELEASE_PADDED_CACHE=1 \
ED2_CSR_DIRECT_BUILD=1 \
ED2_CACHE_COEFF_BOUND=1 \
ED2_CSR_DIRECT_BATCH_ROWS=512 \
./source/QS3ED2 < input.dat
```

The batch setting is optional. If unset, exactly the original **64 rows**
are used. An explicit value is accepted only when direct CSR is
enabled, must consist of decimal digits and lies in the range
**1 through 1024**. Invalid, zero, negative, oversized, nonnumeric
or conflicting settings are rejected. This is an explicitly
user-controlled memory/initialization-time tradeoff; there is
no unvalidated auto-tuning based on nominal OpenMP count.

The native source reports the selected staging width as

    Hamiltonian action: csr-direct-batch-rows=N

and continues to report

    Hamiltonian action: csr-direct-build passes=2, max_staging_bytes=S, peak_cache_payload_upper_bytes=P

The computed upper bound `P` concerns only simultaneous
Hamiltonian CSR+staging array payload, **not process peak RSS**.
For the full matrix of basis dimension `lv`, the stage allocation
has `min(N,lv)` columns. Actual batch rows are smaller for the
last partial batch. The new code bounds its integer indexing
without changing the mathematical H row builder or column order.

All other ED2 modes (including default padded cache, direct flag OFF,
CSR retained and CSR-released-after-conversion) are unaffected by
batch-size control. No other module is permitted to read the padded
arrays after the explicit direct mode or padded-release mode.

## Reproducible numerical proof

The native test `tests/cpc_csr_direct_batch_regression.py`
compares every real and imaginary IEEE-754 bit of **two complex
H*x input/output vectors** against the original 64-row direct
builder for each explicit batch size.

The test covers physically distinct real XYZ, complex DM,
U(1)-preserving XXZ and fully populated exchange tensors, at
N4/N6/N10 (small-system independent dense-H checks at N4/N6),
OpenMP 1/4, batch sizes `1,64,256,512,1024`, under both
`-O2` and historical `-Ofast`. Every batch candidate
must match the default output **bit for bit**, not just to a
numerical tolerance, because there is no mathematical
reassociation between these batch choices.

The test also checks actual native staging byte counts and
total CSR+stage payload, checks byte-identical model input
files, rejects eight invalid/conflicting batch flag cases,
and runs in the existing GNU13/14 comprehensive PR CI,
alongside the separately compiled old-main refactoring
proof and legacy ex1–ex8 checks.

## Exploratory batch performance (NOT PRESTO/CPC timing)

An independent local GNU14.2 `-Ofast -fopenmp` shared-host
pilot measured the full 10-site alternating mixed S1/2,S1
XYZ+DM=0.12 Hamiltonian (D=7776), coefficient bound ON.
For 1 and 4 threads, five fixed batch sizes were each
compiled as **separate source binaries** with no other
mathematical change and compared on identical native input.
All binaries produced the **same full H*x SHA256**.

Six interleaved measured trials per thread/batch combination,
plus one warmup each (60 measured trials total), using
the test-only `system_clock` H-cache constructor timer
and `/usr/bin/time -v` maximum process RSS:

| Row batch | T1 median construction ms | T4 median construction ms | T1 peak RSS MiB | T4 peak RSS MiB | Cache-array payload upper MiB |
|---:|---:|---:|---:|---:|---:|
| 64 (default) | 101.10 | 290.20 | 5.64 | 5.46 | 3.76 |
| 128 | 113.20 | 225.73 | 5.77 | 5.71 | 3.83 |
| 256 | 163.94 | 183.98 | 5.73 | 5.52 | 3.97 |
| 512 | 170.67 | 166.04 | 6.06 | 5.89 | 4.24 |
| 1024 | 231.53 | 156.49 | 6.57 | 6.35 | 4.80 |

This sample suggests larger batches can amortize parallel
region startup overhead under four threads, whereas
single-thread runs may be fastest with smaller batches.
However, the shared machine had substantial contention:
some trials had a large MAD, and the setup cost is
model- and hardware-specific. No stable automatic
selection threshold, statistically controlled speedup,
exclusive HPC allocation, or PRESTO benchmark is claimed.

Because the optional direct mode already favors lower
peak cache memory at the expense of two H-building passes,
the conservative default is **64**. Users with multiple
CPU threads and sufficient memory can choose 256/512/1024
after measuring on their own hardware. Always compare
total solver wall time as well as the initial cache
construction; do not optimize the one-time stage at
the expense of the scientific end-to-end objective.

The original source `-Ofast` padded-versus-CSR last-bit
difference noted in earlier docs is unrelated: **between
two direct CSR batch choices, every output bit must match**,
and that is the requirement checked here.
