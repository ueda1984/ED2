# CPC CPiP author-review gates — ED2 (2026-11-19)

**Status: scientific and editorial preparation; NOT submission-ready.**

The suggested route is **Computer Programs in Physics (CPiP)** of
*Computer Physics Communications*, a software archive article in
the CPC Program Library/Mendeley Data. No editor has approved
this program title, authors, paper type or repository snapshot.
Official scope: https://shop.elsevier.com/journals/computer-physics-communications/0010-4655

## Central contribution and what must not be claimed

ED2 extends the *published* QS³/ED program of Ueda, Yunoki and
Shimokawa (CPC 277, 108369, 2022) from spin-1/2 near-saturation
XXZ and magnon-conserving settings to mixed local spin spaces,
general bilinear anisotropies including oriented DM and symmetric
off-diagonal exchange, Cartesian local fields and static
spin observables, and explicitly defined compatible symmetries.

ED2 source was **already available publicly** in QS-Cube/ED2
before this article. This would be its **first program-description
journal paper**, not a claim that the code or mixed-spin ED
was never publicly available. It does not introduce a new
Lanczos algorithm, a general non-Abelian group interface,
a dynamic structure-factor equivalent to all original
QS³ paths, a new magnetic phase, or an installed library ABI.

The recent **generic cache engineering** (ordered complex CSR,
certified row capacity, optional padding deallocation, two-pass
construction) is a useful *implementation feature* with exact
numerical testing. CSR and two-pass assembly are standard
storage/building strategies; **do not claim their invention**.
Actual hardware-dependent speedup numbers are not yet
qualified for the journal.

## Requirements before submission

| Gate | Evidence already available | Required to close |
| --- | --- | --- |
| **1. Paper scope and approved contributors** | First-program-article draft and predecessor comparison; prior QS³ citation. | Author-approved CPiP thesis, complete authors/order/affiliations, historical contributors and CRediT. |
| **2. Scientific reference and model convention** | Historical 25-point stationary mixed-spin DM reference, Hermiticity/HELLmann-Feynman checks; independent NumPy/SciPy reconstruction; generic mainline H*x exact bit testing in real/complex DM. | Recheck the article's **observable response** on the final frozen executable and reproduce the figures from archived inputs/outputs. The figure checker alone checks transcription. |
| **3. General implementation and reproducibility** | Mainline coefficient capacity, LP64 typed BLAS/LAPACK, CSR and direct two-pass build; GNU13/14 exact native H*x tests and reference dense-H; docs source updated in fork. | Freeze exact code commit and source-archive checksum, confirm all accepted flags/memory claims against that frozen code and reproduce test outputs from the archive. |
| **4. Program package** | Automake executable, source tests, Linux build, standard library quickstart script in `paper/tools/`. | Build/distribute a CPC-compliant program bundle with three approved sample input/output families and a clean independent archive-install check; specify tested Linux/BLAS runtime. |
| **5. Version, licence and ownership** | MIT licence; CITATION.cff exists; published 2022 QS³ predecessor. | Reconcile `configure.ac` **1.0.0**, `CITATION.cff` **0.1.0** and its unverified **2026-01-01** release date; obtain project owner/maintainer confirmation. Do not silently assign a tag or overwrite official software metadata. |
| **6. Upstream and documentation** | Changes merged to fork main; pinned MkDocs strict build; official site currently older. | QS-Cube/ED2 maintainers merge the intended code and docs; verify https://ed2.readthedocs.io/en/latest/ rebuilds from the final official version. |
| **7. Publisher code archive** | Public GitHub and local package plan. | Approved stable archive deposited via CPC Program Library/Mendeley Data; identify real accession/DOI, licences and software availability; no invented DOI. |
| **8. Article PDF and disclosures** | Elsevier `elsarticle` source, separate figure/CSV and memory-table audit scripts, cover-letter draft. | Two-pass LaTeX build; page-wise proof of text, figures/tables and bibliography; approved declarations for funding, conflicts, contributions, data availability and ethics where applicable. |

**Submission must remain gated until 1–8 are resolved.**
This is an editorial/process policy, not a prediction of
acceptance. After review, the authors may decide on a different
article category or journal.

## Performance: not a blanket blocker if claims are scoped

The unmerged PR17 N1000 three-magnon cubic `\Gamma/O_h`
coordinate-cache trial has **48/48 identical audited paired
numerical checks**, but operator-exclusive allocation
attestation is missing. Its medians are exploratory and are
now separated into a development-provenance appendix.
The 2026 generic-mainline cache improvements do not
inherit those PRESTO timing results.

**Do not delay the entire scientific program-paper analysis
solely to obtain new HPC time measurements.** A CPiP
description can make its case through the genuinely
more general *physical model support*, independently
reproducible science, and defensible memory accounting,
while openly declining any hardware-portable speedup claim.
If reviewers or editors require measured performance,
arrange new repeated, controlled trials after a final
source freeze, with raw run logs and exact allocation
metadata. Do not imply their completion now.

If the appendix distracts from the accepted program release,
it can be omitted from the submitted version; its archived
PR17 CSV/audit remains development evidence.

## Current scientific verification boundaries

- **Strong evidence:** published-QS³ fixed-sector and momentum
  continuation on historical source-pinned reference programs;
  independent model construction and full-basis mixed-spin
  DM physics in the earlier manuscript preparation;
  generic H*x bitwise old/new checks and dense-H
  references on latest mainline Fortran.
- **Still needs final-release rerun:** native DM/chirality
  complete end-user response, solver convergence controls,
  reference observables, package builds and snapshots
  on the proposed actual official release code.
- **Not claimed:** an independent all-sector N1000
  three-magnon spectrum, a controlled full-system runtime
  ratio for current mainline, or a validated general
  noncommuting point-group representation in official main.

## Practical build and audit on this paper branch

From the source root:

```bash
python3 -S paper/check_dm_figure.py --self-test
python3 -S paper/check_mainline_cache_claims.py --self-test
python3 -S paper/check_pr17_presto_table.py
python3 -S paper/audit_cpc_submission.py --output cpc-gate-report.json
```

For the independent scientific audit, install NumPy/SciPy and run:

```bash
python3 paper/audit_reference_small_systems.py \
  --csv paper/data/cpc-stationary-dm-reference.csv
```

The source-linked `paper/audit_cpc_submission.py` reports
hard release/approval blockers; `--require-ready` **must
continue to fail** until the final actual approvals, archive
and disclosures are supplied. No automatic tagging,
upstream changes, archival deposit or formal submission
is authorized by a green CI build.
