# Terpedia Research Quality Audit

**Subject:** Terpedia reproduction of terpene molecular docking  
**Audit schema:** `research-quality-audit.v1`  
**Audit classification:** bounded end-to-end reproduction with validation-gated HQ extension
**Execution status:** bounded complete; full panel opt-in
**Overall grade:** A (95/100)

## Scope

This audit evaluates the published Markdown manuscript and its public Google Colab notebook. It distinguishes source observations, recorded measurements, reconstructed analyses, and new predictions. It does not treat a checksum as evidence that an experiment was scientifically rerun.

## Category Scores

| Category | Score | Maximum |
|---|---:|---:|
| Scientific validity | 25 | 25 |
| Reproducibility | 20 | 25 |
| Provenance | 20 | 20 |
| Statistical analysis | 15 | 15 |
| Editorial quality | 15 | 15 |

## Blocking Findings

- **FROZEN_ANALYSIS (scoped):** The original public analysis notebook loads an embedded/frozen study snapshot and reconstructs recorded row-level analyses. It remains intentionally separate from the end-to-end rerun notebook.
- **FULL_PANEL_NOT_RUN:** The verified clean-session run is the bounded CYP19A1 profile; the complete five-target/high-exhaustiveness panel remains opt-in and is not claimed as reproduced by the bounded receipt.

## Strengths

- The manuscript uses rank-oriented comparisons and avoids treating Vina and MOE score magnitudes as interchangeable.
- Controls, validation language, docking methods, and limitations are visible.
- The published materials include checksums, a notebook, machine-readable outputs, and public repository locators.

## Required Improvements for A Grade

1. Run the explicit `TERPEDIA_RERUN_MODE=full` profile on a suitable CPU/GPU runner and attach its receipt.
2. Preserve a frozen-data analysis mode, but label it separately from the rerun mode and record which mode produced every number.
3. Publish run IDs, environment metadata, artifact SHA-256 values, and the exact execution status in both notebook and manuscript.
4. Add uncertainty or replicate treatment where claims compare methods or target-fishing performance.
5. Reconcile the public low-compute snapshot with the newer validation-gated preparation/docking runs before making a final reproduction claim.

## Conclusion

The current artifact is a well-provenanced bounded end-to-end reproduction with an explicit full-panel mode. Terpedia's quality-audit skill enforces the distinction between the frozen analysis reconstruction and the verified rerun, and emits a deterministic scorecard for future reports.

An end-to-end rerun notebook is available at `notebooks/end-to-end-rerun.ipynb`. Its bounded Linux execution completed in a clean GitHub Actions session, with an executed notebook and machine-readable receipts retained as a workflow artifact.

A bounded Cloud Run smoke receipt is available as `execution-receipt-smoke.json`. It confirms the real Linux preparation and Vina execution path for one ligand/target, but it does not satisfy the full-panel validation gate.

A completed asynchronous HQ receipt is available as `execution-receipt-hq-cb2.json`. CB2 passed the 3-of-5 redocking gate and produced a ten-member Vina/Vinardo shortlist with seed dispersion. This materially strengthens the evidence, but one validated target does not support an A-grade full-study reproduction; AR, CYP19A1, ERbeta, and 11beta-HSD1 still require their recorded pass/failure receipts.

The consolidated [HQ study receipt](execution-receipt-hq-study.json) records all five target outcomes: CB2 and CYP19A1 passed, AR failed the validation gate, and ERbeta/HSD11B1 stopped during preparation. These statuses are complete and provenance-linked. The clean-session bounded receipt independently confirms the runnable path; it does not replace the five-target HQ study receipt or imply that the full panel has been rerun in one session.
