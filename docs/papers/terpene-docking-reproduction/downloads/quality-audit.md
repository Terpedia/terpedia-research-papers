# Terpedia Research Quality Audit

**Subject:** Terpedia reproduction of terpene molecular docking  
**Audit schema:** `research-quality-audit.v1`  
**Audit classification:** analysis reconstruction  
**Execution status:** partial  
**Overall grade:** B (85/100)

## Scope

This audit evaluates the published Markdown manuscript and its public Google Colab notebook. It distinguishes source observations, recorded measurements, reconstructed analyses, and new predictions. It does not treat a checksum as evidence that an experiment was scientifically rerun.

## Category Scores

| Category | Score | Maximum |
|---|---:|---:|
| Scientific validity | 25 | 25 |
| Reproducibility | 10 | 25 |
| Provenance | 20 | 20 |
| Statistical analysis | 15 | 15 |
| Editorial quality | 15 | 15 |

## Blocking Findings

- **FROZEN_ANALYSIS:** The original public notebook loads an embedded/frozen study snapshot and reconstructs recorded row-level analyses. It does not itself invoke RDKit, Meeko, PDBFixer, AutoDock Vina, GNINA, or a Terpedia docking job.
- **INCOMPLETE_EXECUTION:** The public Colab rerun path has been published, and asynchronous Cloud Run receipts now cover all five validation outcomes, but the Colab has not yet been executed in a clean session and attached as an executed notebook.

## Strengths

- The manuscript uses rank-oriented comparisons and avoids treating Vina and MOE score magnitudes as interchangeable.
- Controls, validation language, docking methods, and limitations are visible.
- The published materials include checksums, a notebook, machine-readable outputs, and public repository locators.

## Required Improvements for A Grade

1. Make the Colab run end-to-end from a clean session with pinned dependencies, public input structures, explicit seeds, and saved executed outputs.
2. Preserve a frozen-data analysis mode, but label it separately from the rerun mode and record which mode produced every number.
3. Publish run IDs, environment metadata, artifact SHA-256 values, and the exact execution status in both notebook and manuscript.
4. Add uncertainty or replicate treatment where claims compare methods or target-fishing performance.
5. Reconcile the public low-compute snapshot with the newer validation-gated preparation/docking runs before making a final reproduction claim.

## Conclusion

The current artifact is a useful and reasonably well-provenanced analysis reconstruction. It is not yet an A-grade runnable reproduction. Terpedia's quality-audit skill now enforces this distinction and emits a deterministic scorecard for future reports.

An end-to-end rerun notebook has been added at `notebooks/end-to-end-rerun.ipynb`. This audit remains B-grade until that notebook completes in a clean Colab session and its executed output is attached to the report; the presence of executable code is not counted as an executed result.

A bounded Cloud Run smoke receipt is available as `execution-receipt-smoke.json`. It confirms the real Linux preparation and Vina execution path for one ligand/target, but it does not satisfy the full-panel validation gate.

A completed asynchronous HQ receipt is available as `execution-receipt-hq-cb2.json`. CB2 passed the 3-of-5 redocking gate and produced a ten-member Vina/Vinardo shortlist with seed dispersion. This materially strengthens the evidence, but one validated target does not support an A-grade full-study reproduction; AR, CYP19A1, ERbeta, and 11beta-HSD1 still require their recorded pass/failure receipts.

The consolidated [HQ study receipt](execution-receipt-hq-study.json) now records all five target outcomes: CB2 and CYP19A1 passed, AR failed the validation gate, and ERbeta/HSD11B1 stopped during preparation. These statuses are complete and provenance-linked, but the public Colab still needs to be executed from a clean session and the manuscript needs to incorporate the consolidated receipt before the grade can be reconsidered.
