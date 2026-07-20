# Reproduction and Concordance Analysis of Terpene Molecular Docking using AutoDock Vina

## Abstract
This study aimed to reproduce and evaluate the concordance of a published work on terpene target fishing and molecular docking [1]. The original study utilized proprietary MOE GBVI/WSA dG for docking 23 terpenes across five selected protein targets. For this reproduction, an open-source workflow comprising RDKit [6], Meeko [7], and AutoDock Vina [10] was employed. Target fishing frequencies reported in the original paper for CB2, AR, CYP19A1, ERbeta, and 11beta-HSD1 were 18, 15, 13, 11, and 10, respectively. Recomputation from unique supplementary consensus rows yielded frequencies of 19, 13, 11, 10, and 7. The author-curated target fishing set contained 257 normalized pairs, while direct recomputation yielded 222 pairs, with 208 overlapping, resulting in a precision of 0.9369 and recall of 0.8093 against author labels. Notably, beta-myrcene was absent from raw prediction rows but present in the author's consensus and docking data. Cross-engine docking concordance between Vina and MOE, assessed by Spearman rank correlation (rho) and top-five ligand overlap, varied across targets: 11beta-HSD1 (rho 0.559, 4/5 overlap), AR (rho 0.2752, 2/5 overlap), CB2 (rho 0.292, 4/5 overlap), CYP19A1 (rho 0.4435, 4/5 overlap), and ERbeta (rho 0.4161, 5/5 overlap). The original MOE winner was a known ligand for four targets and bornyl acetate for CYP19A1, whereas Vina consistently placed a known ligand first for all five. Nerolidol was the top MOE terpene for four targets, while alpha-bisabolol was the top Vina terpene for all five. Ligand efficiency analysis for CYP19A1 revealed that the original claim of terpenes exceeding known ligands was not reproduced by Vina; Vina mean ligand efficiencies were 0.642 for monoterpenes, 0.440 for sesquiterpenes, and 0.498 for known ligands. An audit of enrichment factors highlighted an internal inconsistency in the reference paper's reported 0.38 value for 11beta-HSD1 sesquiterpenes at a 100% cutoff, which should be 1.0. This reproduction study provides a reproducible sensitivity analysis of the original findings, demonstrating that computational methodology significantly impacts virtual screening outcomes and emphasizing the value of open-source tools for independent validation.

## Keywords
Terpenes, Molecular Docking, AutoDock Vina, MOE, Reproducibility, Target Fishing, Ligand Efficiency, Cannabinoid Receptor 2, Androgen Receptor, CYP19A1, Estrogen Receptor Beta, 11beta-HSD1

## Introduction
Computational methods, particularly molecular docking and target fishing, play a crucial role in modern drug discovery and lead optimization. These techniques enable the prediction of ligand-protein interactions and potential therapeutic targets for small molecules. However, the reproducibility of computational studies, especially those relying on proprietary software and undocumented protocols, remains a significant challenge. Ensuring the transparency and verifiability of such analyses is paramount for scientific rigor and the advancement of chemical biology.

Terpenes, a diverse class of organic compounds abundant in cannabis and other plants, have garnered increasing interest due to their wide array of pharmacological activities. Understanding their molecular interactions with various biological targets is essential for elucidating their therapeutic potential. A recent study investigated the target fishing and molecular docking of 23 terpenes across a panel of five key protein targets using proprietary MOE software [1].

This reproduction study aims to provide an independent, open-source sensitivity analysis of the findings reported in the reference manuscript by Hegde et al. [1]. Our objective is to recapitulate the study's design using publicly available tools (RDKit, Meeko, and AutoDock Vina) and compare the resulting ligand rankings, ligand efficiencies, and target fishing consensus with the original proprietary MOE GBVI/WSA dG results. This effort distinguishes between design recapitulation, source-table reconstruction, and cross-engine numerical replication, contributing to the broader discussion on reproducibility in computational chemistry.

## Materials And Methods

### Study Design
The reproduction study followed the overall design of the reference manuscript [1], focusing on a panel of 23 terpenes and five selected protein targets: Cannabinoid Receptor 2 (CB2, PDB ID: [6KPF](https://www.rcsb.org/structure/6KPF)), Androgen Receptor (AR, PDB ID: [3V49](https://www.rcsb.org/structure/3V49)), Cytochrome P450 19A1 (CYP19A1, PDB ID: [3EQM](https://www.rcsb.org/structure/3EQM)), Estrogen Receptor Beta (ERbeta, PDB ID: [1L2J](https://www.rcsb.org/structure/1L2J)), and 11-beta-hydroxysteroid dehydrogenase 1 (11beta-HSD1, PDB ID: [1XU7](https://www.rcsb.org/structure/1XU7)). The study included co-crystal ligands, as well as high- and low-affinity ligands from the IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb) [11], as controls for each target.

### Target Fishing
The target fishing analysis in the original study utilized four independent prediction engines: SwissTargetPrediction [2], SEA [3], STarFish [5], and PPB2 [4]. A ligand-target pair was considered a consensus prediction if it was predicted by at least a minimum number of independent engines. Raw prediction rows from the author's supplementary information 1 were used to recompute consensus predictions [12].

### Molecular Docking
**Receptor Preparation:** For each target, the protein structure was retrieved from the RCSB Protein Data Bank [6]. Receptors were prepared using Meeko version 0.7.1 [7]. PDBFixer version 1.12.0 was used to add missing atoms and hydrogens at pH 7.4 [8]. Heterogens were excluded from PDBFixer processing and restored unchanged. Local OpenMM minimization with restraints was applied for relaxation [9]. A specific deviation occurred for ERbeta (PDB ID: 1L2J): PDBFixer output created a residue valence conflict in Meeko 0.7.1, so observed chain-A coordinates were passed directly to Meeko, which excluded incomplete residues that could not match a complete chemical template. The reference ligand from the co-crystal structure was removed, and a padding of 5 Å was applied to define the docking box.

**Ligand Preparation:** The 23 terpenes and control ligands were prepared using RDKit ETKDGv3 for conformer generation and Meeko for atom typing and charge assignment [6,7]. PubChem compound identities and structures were used as source records [13].

**Docking Engine:** AutoDock Vina version 1.2.7 with the Vina scoring function was used for all docking calculations [10]. The exhaustiveness parameter was set to 4, and 3 poses were requested per ligand. A consistent random seed (20260717) was used for all runs.

**Computational Execution:** A total of 130 docking jobs were executed. The computations were performed on Google Cloud Run Jobs using the `terpedia-cheminformatics` tool version 1.1.0 (container digest `sha256:95e3067ff1592f17200b72852e9f5e92ef9f1a782f1a960259f5c9b096c1eec6`). This execution strategy involved preparing each receptor once and docking its 26-ligand batch serially on one CPU, with the service scaling to zero between runs. Cloud Run's pay-per-use and scale-to-zero behavior is documented by Google Cloud [14]. This represents a low-cost rank and enrichment robustness check.

### Data Analysis
The primary analysis focused on comparing the independent AutoDock Vina results with the reference MOE GBVI/WSA dG results. This comparison was based solely on within-target ligand ranks and top-set overlap, as score magnitudes are not numerically comparable across different docking engines. Key metrics included:
- **Spearman Rank Correlation (rho):** To assess the correlation of ligand rankings between Vina and MOE for each target.
- **Top-N Overlap:** The number of shared ligands within the top 5 ranked compounds for each target.
- **Ligand Efficiency:** Defined as the absolute Vina score magnitude divided by the number of heavy atoms in the ligand.
- **Enrichment Factor (EF):** Calculated to evaluate the ability of the docking protocol to prioritize known active ligands and TF-positive terpenes.

## Results

### Target Fishing Reproducibility
The target fishing analysis revealed material discrepancies between the reported frequencies in the reference paper and those recomputed from the unique supplementary consensus rows [1,12]. The reported frequencies for CB2, AR, CYP19A1, ERbeta, and 11beta-HSD1 were 18, 15, 13, 11, and 10, respectively. In contrast, the recomputed unique supplementary consensus rows yielded frequencies of 19, 13, 11, 10, and 7 for the same targets.

The author-curated consensus set contained 257 normalized ligand-target pairs. Direct recomputation of consensus from the raw prediction rows resulted in 222 pairs. Comparing these two sets, 208 pairs overlapped, yielding a precision of 0.9369 and a recall of 0.8093 against the author's labels. A notable observation is the absence of beta-myrcene from the raw prediction rows, despite its presence in the author's consensus and docking data [12].

### Cross-Engine Docking Concordance
A comparison of the AutoDock Vina docking results with the reference MOE GBVI/WSA dG results, based on within-target ligand ranks and top-five overlap, is summarized in Table 1 [1,10].

**Table 1: Cross-Engine Docking Concordance (Vina vs. MOE)**

| Target       | Spearman Rho | Top-5 Overlap |
|--------------|--------------|---------------|
| 11beta-HSD1  | 0.559        | 4/5           |
| AR           | 0.2752       | 2/5           |
| CB2          | 0.292        | 4/5           |
| CYP19A1      | 0.4435       | 4/5           |
| ERbeta       | 0.4161       | 5/5           |

For the top-ranked ligands, the original MOE analysis identified a known ligand as the top binder for four of the five targets (11beta-HSD1, AR, CB2, ERbeta), with bornyl acetate being the top ligand for CYP19A1. In contrast, the Vina reproduction consistently placed a known ligand first for all five targets. Among the terpenes, nerolidol was identified as the top MOE terpene for four targets, while alpha-bisabolol emerged as the top Vina terpene for all five targets.

### Ligand Efficiency Analysis
The ligand efficiency analysis for CYP19A1 demonstrated a notable difference from the original study's claims [1]. The original paper suggested that every terpene exceeded every known ligand in terms of ligand efficiency for CYP19A1. However, the Vina reproduction did not support this claim. The mean ligand efficiencies calculated from the Vina scores for CYP19A1 were 0.642 for monoterpenes, 0.440 for sesquiterpenes, and 0.498 for known ligands [10,16].

### Enrichment Factor Audit
A critical audit of the enrichment factor (EF) definition revealed an internal inconsistency in the reference paper [1]. Under the stated EF definition, the EF at a 100 percent cutoff must theoretically equal 1.0. However, the paper reported an EF of 0.38 for 11beta-HSD1 sesquiterpenes at a 100% cutoff, which is inconsistent with the definition and therefore not considered valid in this reproduction [16].

## Discussion
This reproduction study provides a sensitivity analysis of the target fishing and molecular docking results presented in the reference manuscript [1], using an open-source computational workflow [6,7,10]. The findings highlight both areas of concordance and significant differences when comparing results from proprietary MOE software with those from RDKit/Meeko/AutoDock Vina.

The observed Spearman rank correlations between the two docking engines range from low (AR: 0.2752, CB2: 0.292) to modest (ERbeta: 0.4161, CYP19A1: 0.4435) and moderate (11beta-HSD1: 0.559). This suggests that while there is some agreement in the overall ranking of ligands, particularly for 11beta-HSD1, the specific order and relative potency predictions can vary considerably depending on the docking engine and its underlying scoring function. The top-five overlap analysis serves as an exploratory measure of cross-engine agreement, indicating the extent to which the most promising candidates are consistently identified. However, it is important to note that this overlap does not establish experimental binding affinity, predictive accuracy, or biological validity.

The differences in the top-ranked ligands and the ligand efficiency analysis for CYP19A1 underscore the methodological sensitivity of virtual screening campaigns. The original claim regarding terpene ligand efficiency for CYP19A1 was not reproduced by the Vina workflow, suggesting that the choice of docking software and scoring function can profoundly influence the interpretation of results, particularly for novel chemical spaces like terpenes.

The target fishing reproducibility section revealed material discrepancies in target frequencies and consensus predictions. The absence of beta-myrcene from the raw prediction rows, despite its inclusion in the author's consensus and docking, points to potential data handling or reporting inconsistencies in the original supplementary materials. These findings emphasize the importance of transparent and fully accessible raw data for complete reproducibility.

## Reference Figure and Table Reconstruction
The public analysis Colab executes an analysis reconstruction from the recorded study snapshot. A separate [end-to-end rerun Colab](../notebooks/end-to-end-rerun.ipynb) now installs the pinned chemistry environment, runs receptor preparation and validation-gated docking, and writes a machine-readable receipt; its execution must be completed and its receipt attached before the study can claim a complete rerun. No source figures were supplied in the reference materials available to this study, so no figure is claimed as reproduced.

| Source display | Terpedia display | Status | Artifact or blocking reason |
|---|---|---|---|
| Table 1 | Table 1 | reproduced | `results.xlsx`, analysis cells in `reproducibility.ipynb` |
| Reference figures | None supplied | blocked | Source figure files were not present in the supplied reference package |

This distinction is part of the reproducibility record and is not evidence of successful end-to-end replication.

## Limitations
This reproduction study has several limitations. The primary limitation is that numerical values of docking scores from the proprietary MOE GBVI/WSA dG method are not directly comparable with those from AutoDock Vina. Therefore, comparisons were restricted to ligand ranks and top-set overlaps. While row-level records from all four target-fishing tools were available and Terpedia computed leave-one-tool-out ablation, the target-fishing analysis was limited by the disagreement between author-curated and directly recomputed consensus, uneven tool coverage, and unresolved source-label reconciliation. The absence of beta-myrcene from the raw prediction rows, despite its presence in the author's consensus and docking, represents a data gap. Furthermore, the receptor preparation for ERbeta involved a deviation from the standard protocol due to a residue valence conflict in Meeko 0.7.1, which may impact the generalizability of results for this specific target. The top-five overlap is an exploratory measure of cross-engine agreement and does not establish binding, affinity, predictive accuracy, or experimental validity.

## Conclusion
This study implemented an open-source RDKit/Meeko/AutoDock Vina workflow to perform a reproducible sensitivity analysis of the target fishing and molecular docking results presented in the reference manuscript. The public notebook is an analysis reconstruction of recorded results, not an end-to-end rerun of the docking experiment. While a moderate rank concordance was observed for 11beta-HSD1 (Spearman rho 0.559) and modest concordance for ERbeta (Spearman rho 0.4161) and CYP19A1 (Spearman rho 0.4435), lower correlations were found for AR (Spearman rho 0.2752) and CB2 (Spearman rho 0.292). Significant differences were identified in top-ranked ligands, ligand efficiency for CYP19A1, and target fishing consensus. These findings collectively underscore the critical impact of computational methodology on virtual screening outcomes and highlight the value of open-source tools for conducting independent sensitivity analyses and fostering greater transparency in computational drug discovery.

## Data Availability
The data and computational artifacts supporting this study are available in the public [Terpedia research-papers repository](https://github.com/Terpedia/terpedia-research-papers/tree/main/docs/papers/terpene-docking-reproduction) [15], including the [reproducibility notebook](https://github.com/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/notebooks/reproducibility.ipynb), [results workbook](https://github.com/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/downloads/results.xlsx), and [target-fishing records](https://github.com/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/downloads/reference-target-fishing.json). The author-supplied reference manuscript and supplements are preserved through the [source manifest](https://github.com/Terpedia/chat-terpedia-backend/blob/main/research/reference-paper/supplementary-manifest.json) and [author-supplied Drive folder](https://drive.google.com/drive/folders/1ooJl0fn4EGIy7AHjs52qHYSG9-nqRBs-); these source links may require authorization. The specific study bundle has a SHA-256 checksum of `cc88f76d5c344e64dc85f5210adca14bb9a6dac53b6d7421a790a38730d14020`. The notebook is also available in [Google Colab](https://colab.research.google.com/github/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/notebooks/reproducibility.ipynb). The notebook content SHA-256 is `d863b529d17d85034e6af6f9a430526a23127d80fbe9600fa9f89d1c6786d52f`, and the underlying research snapshot SHA-256 is `a0f2569040758f879017e98239620e7a280b83b1830c73f5d3a74545f5b7e325`. The Terpedia chemistry service version used for this analysis was 1.2.0. Docking calculations were performed on Google Cloud Run Jobs using the `terpedia-cheminformatics` tool version 1.1.0, with a container digest of `sha256:95e3067ff1592f17200b72852e9f5e92ef9f1a782f1a960259f5c9b096c1eec6`.

## AI Disclosure
This manuscript was generated and edited by the Terpedia Editor Agent, an AI system, with engineering assistance from Codex.

## References
[1] Hegde M, Ahmadi A, McShan D, Trapp S, Baudry J. *Target Fishing and Molecular Docking of Terpenes: A Systematic Evaluation of Method Concordance and Ligand Efficiency*. Author-supplied manuscript, version received 2026, with supplements listed in the [Terpedia source manifest](https://github.com/Terpedia/chat-terpedia-backend/blob/main/research/reference-paper/supplementary-manifest.json).
[2] Daina A, Michielin O, Zoete V. SwissTargetPrediction: updated data and new features for efficient prediction of protein targets of small molecules. *Nucleic Acids Research*. 2019;47(W1):W357-W364. doi:[10.1093/nar/gkz382](https://doi.org/10.1093/nar/gkz382).
[3] Keiser MJ, Roth BL, Armbruster BN, Ernsberger P, Irwin JJ, Shoichet BK. Predicting new molecular targets for known drugs. *Nature*. 2009;462:175-181. doi:[10.1038/nature07839](https://doi.org/10.1038/nature07839).
[4] Awale M, Reymond JL. Polypharmacology Browser PPB2: Target Prediction Combining Nearest Neighbors with Machine Learning. *Journal of Chemical Information and Modeling*. 2019;59(1):10-17. doi:[10.1021/acs.jcim.8b00524](https://doi.org/10.1021/acs.jcim.8b00524).
[5] Cockroft NT, Cheng X, Fuchs JR. STarFish: A Stacked Ensemble Target Fishing Approach and its Application to Natural Products. *Journal of Chemical Information and Modeling*. 2019;59(11):4906-4920. doi:[10.1021/acs.jcim.9b00489](https://doi.org/10.1021/acs.jcim.9b00489).
[6] Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids Research*. 2000;28(1):235-242. doi:[10.1093/nar/28.1.235](https://doi.org/10.1093/nar/28.1.235); Landrum G. RDKit: Open-source cheminformatics. Version 2025.09.5. doi:[10.5281/zenodo.591637](https://doi.org/10.5281/zenodo.591637).
[7] Santos-Martins D, He Y, Eberhardt J, et al. Meeko: Molecule Parametrization and Software Interoperability for Docking and Beyond. *Journal of Chemical Information and Modeling*. 2025;65(24):13045-13050. doi:[10.1021/acs.jcim.5c02271](https://doi.org/10.1021/acs.jcim.5c02271).
[8] Eastman P, Pande VS. PDBFixer: a Python script for fixing problems in PDB files. OpenMM project, version 1.12. [Source repository](https://github.com/openmm/pdbfixer).
[9] Eastman P, Galvelis R, Peláez RP, et al. OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials. *PLoS Computational Biology*. 2024;20(1):e1011460. doi:[10.1371/journal.pcbi.1011460](https://doi.org/10.1371/journal.pcbi.1011460).
[10] Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *Journal of Chemical Information and Modeling*. 2021;61(8):3891-3898. doi:[10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203).
[11] Harding SD, Faccenda E, Southan C, et al. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. *Nucleic Acids Research*. 2024;52(D1):D1327-D1334. doi:[10.1093/nar/gkad944](https://doi.org/10.1093/nar/gkad944).
[12] Terpedia. Author-supplied supplementary-information snapshot and derived target-fishing records. Version retrieved 2026-07-19. [Manifest and checksums](https://github.com/Terpedia/chat-terpedia-backend/blob/main/research/reference-paper/supplementary-manifest.json).
[13] Kim S, Chen J, Cheng T, et al. PubChem 2019 update: improved access to chemical data. *Nucleic Acids Research*. 2019;47(D1):D1102-D1109. doi:[10.1093/nar/gky1033](https://doi.org/10.1093/nar/gky1033).
[14] Google Cloud. Cloud Run pricing and billing documentation. Accessed 2026-07-19. [Cloud Run pricing](https://cloud.google.com/run/pricing).
[15] Terpedia. Reproducibility notebook and public analysis artifacts for this study. Version published 2026-07-19. [Public paper repository](https://github.com/Terpedia/terpedia-research-papers/tree/main/docs/papers/terpene-docking-reproduction).
[16] Terpedia. Machine-readable docking observations, consensus rows, and analysis tables. Checksummed study artifacts published with this paper; see [results.xlsx](https://github.com/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/downloads/results.xlsx) and [reference-target-fishing.json](https://github.com/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/downloads/reference-target-fishing.json).
