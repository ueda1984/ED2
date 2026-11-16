# Legacy QS³ versus ED2 at large site counts: experiment

The historical QS-Cube/ED source is pinned at
`a3b09934bba5f966bfb564718528de8778c2b46c`.
The modern ED2 version is the GitHub Actions checkout SHA.
Source is unmodified; both use fixed Ndown, identical XXZ
couplings, the same momentum, and native symmetry-adapted bases.

Default CI (seven source-pinned regressions): 1D periodic rings
L=12,Ndown=2,k=0; L=24,Ndown=3,k=0; L=48,Ndown=2,k=1;
L=96,Ndown=2,k=0; cubic 3×3×3 and 4×4×4 Ndown=2 at Gamma;
cubic 10×10×10 Ndown=1 at R. Shared-runner wall time and RSS
are archived but are **not** an HPC/node-controlled performance claim.

A separate 10×10×10, 1000-site, Ndown=2, Gamma-point stress
with 330-second per-executable timeout is explicitly **non-gating**.
Failure, timeout and asymmetry between executables must be reported
honestly; neither completion nor three-magnon continuity is assumed.

The independent sector dimension is computed by character-weighted
Burnside counting. A translation of order t has C(L/t,Ndown/t)
fixed bitstrings if t divides Ndown, or zero otherwise.
At 10×10×10 Gamma, the expected n=1,2,3 sector dimensions
are 1, 503 and 166167 respectively. The last number is a
**combinatorial prediction, not a validated ED2 run**.
For n=1 the energy is checked against the analytic lattice dispersion.

Original QS³ and modern ED2 use independent random Lanczos
initialization. Energy agreement does not imply matched
iteration histories; observed times are exploratory and not speedups.
A QuSpin comparison must use the same Hamiltonian, conservation
sector, momentum, hardware and accuracy criteria.
