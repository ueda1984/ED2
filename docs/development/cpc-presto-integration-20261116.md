# OpenMP integration and PRESTO measurement checkpoint

The approved OpenMP implementation is incorporated into the native source, not applied at runtime.
Baseline: `719da7035e66ce884e9de00376e4ded5b3e5c459`.
Review source: `e042e90a1d985b2c8320f19d9de2e18d779a273e`.
Native source hashes remain exactly those verified by the OpenMP audit.

This checkpoint freezes the algorithm for the first PRESTO measurements. No further speculative performance changes, reduced precision, relaxed residual tolerance or global thread limit are introduced.

PRESTO stages: (1) read-only queue/environment checks; (2) debug-queue shared four-core smoke, not performance evidence; (3) operator-approved single-host exclusive standard-queue comparison, with physical-core binding, BLAS=1, two warmups and seven measured repeats. Preserve existing jobs, queue policies and fair share. Do not use the closed workq or infer node exclusivity from ncpus alone.

Keep an immutable source archive, compiler/build records, binary and input hashes, full solver diagnostics, individual timing/RSS records and failure records. Compare old/new ED2 with the same full tight cache first; measure the new version's cache fractions separately. External-library timing is a separate, still unfinished experiment.

The submission manuscript remains provisional: shared-host timings are not exclusive-node publication measurements, 16/32-thread scaling has not yet been established, and no PRESTO job has been submitted by this integration workflow.
