# Reproduction and Concordance Analysis of Terpene Molecular Docking using AutoDock Vina

## Abstract
This study aimed to reproduce and evaluate the concordance of a published work on terpene target fishing and molecular docking. The original study utilized proprietary MOE GBVI/WSA dG for docking 23 terpenes across five selected protein targets. For this reproduction, an open-source workflow comprising RDKit, Meeko, and AutoDock Vina was employed. Target fishing frequencies reported in the original paper for CB2, AR, CYP19A1, ERbeta, and 11beta-HSD1 were 18, 15, 13, 11, and 10, respectively. Recomputation from unique supplementary consensus rows yielded frequencies of 19, 13, 11, 10, and 7. The author-curated target fishing set contained 257 normalized pairs, while direct recomputation yielded 222 pairs, with 208 overlapping, resulting in a precision of 0.9369 and recall of 0.8093 against author labels. Notably, beta-myrcene was absent from raw prediction rows but present in the author's consensus and docking data. Cross-engine docking concordance between Vina and MOE, assessed by Spearman rank correlation (rho) and top-five ligand overlap, varied across targets: 11beta-HSD1 (rho 0.559, 4/5 overlap), AR (rho 0.2752, 2/5 overlap), CB2 (rho 0.292, 4/5 overlap), CYP19A1 (rho 0.4435, 4/5 overlap), and ERbeta (rho 0.4161, 5/5 overlap). The original MOE winner was a known ligand for four targets and bornyl acetate for CYP19A1, whereas Vina consistently placed a known ligand first for all five. Nerolidol was the top MOE terpene for four targets, while alpha-bisabolol was the top Vina terpene for all five. Ligand efficiency analysis for CYP19A1 revealed that the original claim of terpenes exceeding known ligands was not reproduced by Vina; Vina mean ligand efficiencies were 0.642 for monoterpenes, 0.440 for sesquiterpenes, and 0.498 for known ligands. An audit of enrichment factors highlighted an internal inconsistency in the reference paper's reported 0.38 value for 11beta-HSD1 sesquiterpenes at a 100% cutoff, which should be 1.0. This reproduction study provides a reproducible sensitivity analysis of the original findings, demonstrating that computational methodology significantly impacts virtual screening outcomes and emphasizing the value of open-source tools for independent validation.

## Keywords
Terpenes, Molecular Docking, AutoDock Vina, MOE, Reproducibility, Target Fishing, Ligand Efficiency, Cannabinoid Receptor 2, Androgen Receptor, CYP19A1, Estrogen Receptor Beta, 11beta-HSD1

## Introduction
Computational methods, particularly molecular docking and target fishing, play a crucial role in modern drug discovery and lead optimization. These techniques enable the prediction of ligand-protein interactions and potential therapeutic targets for small molecules. However, the reproducibility of computational studies, especially those relying on proprietary software and undocumented protocols, remains a significant challenge. Ensuring the transparency and verifiability of such analyses is paramount for scientific rigor and the advancement of chemical biology.

Terpenes, a diverse class of organic compounds abundant in cannabis and other plants, have garnered increasing interest due to their wide array of pharmacological activities. Understanding their molecular interactions with various biological targets is essential for elucidating their therapeutic potential. A recent study investigated the target fishing and molecular docking of 23 terpenes across a panel of five key protein targets using proprietary MOE software.

This reproduction study aims to provide an independent, open-source sensitivity analysis of the findings reported in the reference manuscript by Hegde et al. [1]. Our objective is to recapitulate the study's design using publicly available tools (RDKit, Meeko, and AutoDock Vina) and compare the resulting ligand rankings, ligand efficiencies, and target fishing consensus with the original proprietary MOE GBVI/WSA dG results. This effort distinguishes between design recapitulation, source-table reconstruction, and cross-engine numerical replication, contributing to the broader discussion on reproducibility in computational chemistry.

## Materials And Methods

### Study Design
The reproduction study followed the overall design of the reference manuscript [1], focusing on a panel of 23 terpenes and five selected protein targets: Cannabinoid Receptor 2 (CB2, PDB ID: 6KPF), Androgen Receptor (AR, PDB ID: 3V49), Cytochrome P450 19A1 (CYP19A1, PDB ID: 3EQM), Estrogen Receptor Beta (ERbeta, PDB ID: 1L2J), and 11-beta-hydroxysteroid dehydrogenase 1 (11beta-HSD1, PDB ID: 1XU7). The study included co-crystal ligands, as well as high- and low-affinity ligands from GtoPdb, as controls for each target.

### Target Fishing
The target fishing analysis in the original study utilized four independent prediction engines: SwissTargetPrediction, SEA, STarFish, and PPB2. A ligand-target pair was considered a consensus prediction if it was predicted by at least a minimum number of independent engines. Raw prediction rows from the author's supplementary information 1 were used to recompute consensus predictions.

### Molecular Docking
**Receptor Preparation:** For each target, the protein structure was retrieved from the RCSB Protein Data Bank. Receptors were prepared using Meeko version 0.7.1. PDBFixer version 1.12.0 was used to add missing atoms and hydrogens at pH 7.4. Heterogens were excluded from PDBFixer processing and restored unchanged. Local OpenMM minimization with restraints was applied for relaxation. A specific deviation occurred for ERbeta (PDB ID: 1L2J): PDBFixer output created a residue valence conflict in Meeko 0.7.1, so observed chain-A coordinates were passed directly to Meeko, which excluded incomplete residues that could not match a complete chemical template. The reference ligand from the co-crystal structure was removed, and a padding of 5 Å was applied to define the docking box.

**Ligand Preparation:** The 23 terpenes and control ligands were prepared using RDKit ETKDGv3 for conformer generation and Meeko for atom typing and charge assignment.

**Docking Engine:** AutoDock Vina version 1.2.7 with the Vina scoring function was used for all docking calculations. The exhaustiveness parameter was set to 4, and 3 poses were requested per ligand. A consistent random seed (20260717) was used for all runs.

**Computational Execution:** A total of 130 docking jobs were executed. The computations were performed on Google Cloud Run Jobs using the `terpedia-cheminformatics` tool version 1.1.0 (container digest `sha256:95e3067ff1592f17200b72852e9f5e92ef9f1a782f1a960259f5c9b096c1eec6`). This execution strategy involved preparing each receptor once and docking its 26-ligand batch serially on one CPU, with the service scaling to zero between runs. This represents a low-cost rank and enrichment robustness check.

### Data Analysis
The primary analysis focused on comparing the independent AutoDock Vina results with the reference MOE GBVI/WSA dG results. This comparison was based solely on within-target ligand ranks and top-set overlap, as score magnitudes are not numerically comparable across different docking engines. Key metrics included:
- **Spearman Rank Correlation (rho):** To assess the correlation of ligand rankings between Vina and MOE for each target.
- **Top-N Overlap:** The number of shared ligands within the top 5 ranked compounds for each target.
- **Ligand Efficiency:** Defined as the absolute Vina score magnitude divided by the number of heavy atoms in the ligand.
- **Enrichment Factor (EF):** Calculated to evaluate the ability of the docking protocol to prioritize known active ligands and TF-positive terpenes.

## Results

### Target Fishing Reproducibility
The target fishing analysis revealed material discrepancies between the reported frequencies in the reference paper and those recomputed from the unique supplementary consensus rows. The reported frequencies for CB2, AR, CYP19A1, ERbeta, and 11beta-HSD1 were 18, 15, 13, 11, and 10, respectively. In contrast, the recomputed unique supplementary consensus rows yielded frequencies of 19, 13, 11, 10, and 7 for the same targets.

The author-curated consensus set contained 257 normalized ligand-target pairs. Direct recomputation of consensus from the raw prediction rows resulted in 222 pairs. Comparing these two sets, 208 pairs overlapped, yielding a precision of 0.9369 and a recall of 0.8093 against the author's labels. A notable observation is the absence of beta-myrcene from the raw prediction rows, despite its presence in the author's consensus and docking data.

### Cross-Engine Docking Concordance
A comparison of the AutoDock Vina docking results with the reference MOE GBVI/WSA dG results, based on within-target ligand ranks and top-five overlap, is summarized in Table 1.

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
The ligand efficiency analysis for CYP19A1 demonstrated a notable difference from the original study's claims. The original paper suggested that every terpene exceeded every known ligand in terms of ligand efficiency for CYP19A1. However, the Vina reproduction did not support this claim. The mean ligand efficiencies calculated from the Vina scores for CYP19A1 were 0.642 for monoterpenes, 0.440 for sesquiterpenes, and 0.498 for known ligands.

### Enrichment Factor Audit
A critical audit of the enrichment factor (EF) definition revealed an internal inconsistency in the reference paper. Under the stated EF definition, the EF at a 100 percent cutoff must theoretically equal 1.0. However, the paper reported an EF of 0.38 for 11beta-HSD1 sesquiterpenes at a 100% cutoff, which is inconsistent with the definition and therefore not considered valid in this reproduction.

## Discussion
This reproduction study provides a sensitivity analysis of the target fishing and molecular docking results presented in the reference manuscript [1], using an open-source computational workflow. The findings highlight both areas of concordance and significant differences when comparing results from proprietary MOE software with those from RDKit/Meeko/AutoDock Vina.

The observed Spearman rank correlations between the two docking engines range from low (AR: 0.2752, CB2: 0.292) to modest (ERbeta: 0.4161, CYP19A1: 0.4435) and moderate (11beta-HSD1: 0.559). This suggests that while there is some agreement in the overall ranking of ligands, particularly for 11beta-HSD1, the specific order and relative potency predictions can vary considerably depending on the docking engine and its underlying scoring function. The top-five overlap analysis serves as an exploratory measure of cross-engine agreement, indicating the extent to which the most promising candidates are consistently identified. However, it is important to note that this overlap does not establish experimental binding affinity, predictive accuracy, or biological validity.

The differences in the top-ranked ligands and the ligand efficiency analysis for CYP19A1 underscore the methodological sensitivity of virtual screening campaigns. The original claim regarding terpene ligand efficiency for CYP19A1 was not reproduced by the Vina workflow, suggesting that the choice of docking software and scoring function can profoundly influence the interpretation of results, particularly for novel chemical spaces like terpenes.

The target fishing reproducibility section revealed material discrepancies in target frequencies and consensus predictions. The absence of beta-myrcene from the raw prediction rows, despite its inclusion in the author's consensus and docking, points to potential data handling or reporting inconsistencies in the original supplementary materials. These findings emphasize the importance of transparent and fully accessible raw data for complete reproducibility.

## Limitations
This reproduction study has several limitations. The primary limitation is that numerical values of docking scores from the proprietary MOE GBVI/WSA dG method are not directly comparable with those from AutoDock Vina. Therefore, comparisons were restricted to ligand ranks and top-set overlaps. While row-level records from all four target-fishing tools were available and Terpedia computed leave-one-tool-out ablation, the target-fishing analysis was limited by the disagreement between author-curated and directly recomputed consensus, uneven tool coverage, and unresolved source-label reconciliation. The absence of beta-myrcene from the raw prediction rows, despite its presence in the author's consensus and docking, represents a data gap. Furthermore, the receptor preparation for ERbeta involved a deviation from the standard protocol due to a residue valence conflict in Meeko 0.7.1, which may impact the generalizability of results for this specific target. The top-five overlap is an exploratory measure of cross-engine agreement and does not establish binding, affinity, predictive accuracy, or experimental validity.

## Conclusion
This reproduction study successfully implemented an open-source RDKit/Meeko/AutoDock Vina workflow to re-evaluate the target fishing and molecular docking of terpenes as presented in the reference manuscript. The workflow provides a reproducible sensitivity analysis of the original proprietary MOE GBVI/WSA dG results. While a moderate rank concordance was observed for 11beta-HSD1 (Spearman rho 0.559) and modest concordance for ERbeta (Spearman rho 0.4161) and CYP19A1 (Spearman rho 0.4435), lower correlations were found for AR (Spearman rho 0.2752) and CB2 (Spearman rho 0.292). Significant differences were identified in top-ranked ligands, ligand efficiency for CYP19A1, and target fishing consensus. These findings collectively underscore the critical impact of computational methodology on virtual screening outcomes and highlight the value of open-source tools for conducting independent sensitivity analyses and fostering greater transparency in computational drug discovery.

## Data Availability
The data and computational artifacts supporting this study are available in the Terpedia research repository: `https://github.com/Terpedia/terpedia-reproduce-the-research-in-the-attached-reference-manusc-be6fa4a4`. The specific study bundle has a SHA-256 checksum of `cc88f76d5c344e64dc85f5210adca14bb9a6dac53b6d7421a790a38730d14020`. The complete reproducibility notebook, including all code and analysis, is publicly available at: `https://colab.research.google.com/github/Terpedia/terpedia-research-papers/blob/main/docs/papers/terpene-docking-reproduction/notebooks/reproducibility.ipynb`. The notebook content SHA-256 is `d863b529d17d85034e6af6f9a430526a23127d80fbe9600fa9f89d1c6786d52f`, and the underlying research snapshot SHA-256 is `a0f2569040758f879017e98239620e7a280b83b1830c73f5d3a74545f5b7e325`. The Terpedia chemistry service version used for this analysis was 1.2.0. Docking calculations were performed on Google Cloud Run Jobs using the `terpedia-cheminformatics` tool version 1.1.0, with a container digest of `sha256:95e3067ff1592f17200b72852e9f5e92ef9f1a782f1a960259f5c9b096c1eec6`.

## AI Disclosure
This manuscript was generated and edited by the Terpedia Editor Agent, an AI system, with engineering assistance from Codex.

## References
[1] Hegde M, Ahmadi A, McShan D, Trapp S, Baudry J. Supplied manuscript, version received 2026.