# Terpedia Research Quality Audit

**Subject:** Terpedia reproduction of terpene molecular docking  
**Audit schema:** `research-quality-audit.v1`  
**Audit classification:** analysis reconstruction  
**Execution status:** partial  
**Overall grade:** C (80/100)

## Scope

This audit evaluates the published Markdown manuscript and its public Google Colab notebook. It distinguishes source observations, recorded measurements, reconstructed analyses, and new predictions. It does not treat a checksum as evidence that an experiment was scientifically rerun.

## Category Scores

| Category | Score | Maximum |
|---|---:|---:|
| Scientific validity | 25 | 25 |
| Reproducibility | 10 | 25 |
| Provenance | 20 | 20 |
| Statistical analysis | 10 | 15 |
| Editorial quality | 15 | 15 |

## Blocking Findings

- **FROZEN_ANALYSIS:** The notebook loads an embedded/frozen study snapshot and reconstructs recorded row-level analyses. It does not invoke RDKit, Meeko, PDBFixer, AutoDock Vina, GNINA, or a Terpedia docking job to rerun the experiment.
- **INCOMPLETE_EXECUTION:** The public artifact is therefore a partial analysis reconstruction, not an end-to-end reproduction. The paper must preserve that distinction in its title, abstract, results, and limitations.
- **UNCERTAINTY:** The published analysis does not expose uncertainty, replicate, bootstrap, or confidence-interval treatment for its comparative claims.

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

An end-to-end rerun notebook has been added at `notebooks/end-to-end-rerun.ipynb`. This audit remains C-grade until that notebook completes in a clean Colab session and its JSON receipt is attached to the report; the presence of executable code is not counted as an executed result.
