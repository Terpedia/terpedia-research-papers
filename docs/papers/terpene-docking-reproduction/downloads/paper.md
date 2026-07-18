# Reproduction and Concordance Analysis of Terpene Molecular Docking using AutoDock Vina

## Abstract
**Background:** Terpenes, a diverse class of natural products, are increasingly recognized for their potential therapeutic properties, often mediated through interactions with various biological targets. Computational methods, such as molecular docking, are crucial for identifying and characterizing these interactions. A previous study utilized proprietary MOE software for target fishing and molecular docking of 23 terpenes against a panel of five targets. This work aims to reproduce and evaluate the rank concordance of the docking component of that study using a fully open-source and reproducible workflow.

**Methods:** We developed and executed an open-source molecular docking protocol using RDKit, Meeko, and AutoDock Vina. This protocol was applied to the same set of 23 terpenes and relevant control ligands against the five specified protein targets: Cannabinoid Receptor 2 (CB2/6KPF), Androgen Receptor (AR/3V49), Cytochrome P450 19A1 (CYP19A1/3EQM), Estrogen Receptor Beta (ERbeta/1L2J), and 11-beta-Hydroxysteroid Dehydrogenase 1 (11beta-HSD1/1XU7). Ligand efficiency was calculated as the absolute Vina docking score magnitude divided by the number of heavy atoms. Rank concordance between the newly computed Vina docking results and the reference MOE results was assessed using Spearman rank correlation coefficients and top-5 ligand overlap (Jaccard index) for each target. Due to the unavailability of row-level target-fishing outputs from the original study, target-fishing consensus, TF+/TF- enrichment, and tool ablation analyses could not be recomputed.

**Results:** A total of 130 Vina docking observations were successfully generated across the five targets. Comparisons were made for 129 shared ligand-target pairs between the Vina and reference MOE docking results. Spearman rank correlations for the five targets were: 11beta-HSD1 (0.559, moderate), AR (0.2752, low), CB2 (0.292, low), CYP19A1 (0.4435, modest), and ERbeta (0.4161, modest). The top-5 ligand overlaps (Jaccard index) were: 11beta-HSD1 (0.6667, 4/6 ligands), AR (0.25, 2/8 ligands), CB2 (0.6667, 4/6 ligands), CYP19A1 (0.6667, 4/6 ligands), and ERbeta (1.0, 5/5 ligands). The top-5 ligand overlap is an exploratory measure of cross-engine agreement and does not establish binding, affinity, predictive accuracy, or experimental validity.

**Conclusion:** The open-source RDKit/Meeko/AutoDock Vina protocol demonstrated varying degrees of rank concordance with the proprietary MOE method across the five targets. These findings highlight the utility of open-source tools for providing a reproducible sensitivity analysis of computational drug discovery efforts, while also underscoring the impact of docking engine choice and receptor preparation on predicted binding profiles.

## Keywords
Terpenes, molecular docking, AutoDock Vina, MOE, reproduction study, rank correlation, ligand efficiency, cannabinoid receptor 2, androgen receptor, cytochrome P450 19A1, estrogen receptor beta, 11-beta-hydroxysteroid dehydrogenase 1.

## Introduction
The cannabis plant produces a rich array of secondary metabolites, among which terpenes are prominent for their diverse biological activities and contributions to the plant's aroma and flavor [1]. Beyond their sensory roles, terpenes have garnered significant scientific interest due to their potential therapeutic applications, including anti-inflammatory, anxiolytic, and anti-cancer effects [2,3]. Understanding the molecular mechanisms underlying these effects often involves identifying and characterizing the interactions between terpenes and specific biological targets, such as G-protein coupled receptors, nuclear receptors, and enzymes.

Computational approaches, particularly molecular docking, offer a cost-effective and efficient means to predict these ligand-target interactions. Molecular docking simulates the binding of small molecules (ligands) into the active site of a macromolecular target, providing insights into potential binding modes and relative binding affinities. A previous study employed proprietary MOE software for a comprehensive analysis, including target fishing and molecular docking of 23 terpenes against a panel of five key biological targets. Given the importance of reproducibility and accessibility in scientific research, this work aims to re-evaluate the docking component of that study using a fully open-source and transparent computational workflow. Our objective is to assess the rank concordance between the open-source AutoDock Vina results and the proprietary MOE results, thereby providing a sensitivity analysis of the original findings under different computational methodologies.

## Materials and Methods

### Study Design
This reproduction study focused on the molecular docking component of the original research. The primary goal was to compare the ranking of ligands by an open-source docking engine (AutoDock Vina) against the rankings reported from the proprietary MOE software for a shared set of terpenes and control ligands across five protein targets. Due to the unavailability of row-level outputs from the four target-fishing engines used in the original study, it was not possible to recompute target-fishing consensus, TF+/TF- enrichment, or tool ablation analyses.

### Ligand and Target Selection
The study utilized the same panel of 23 terpenes and a set of known control ligands for each target as in the original study. The five protein targets investigated were:
*   Cannabinoid Receptor 2 (CB2, PDB ID: 6KPF, UniProt ID: P34972)
*   Androgen Receptor (AR, PDB ID: 3V49, UniProt ID: P10275)
*   Cytochrome P450 19A1 (CYP19A1, PDB ID: 3EQM, UniProt ID: P11511)
*   Estrogen Receptor Beta (ERbeta, PDB ID: 1L2J, UniProt ID: Q92731)
*   11-beta-Hydroxysteroid Dehydrogenase 1 (11beta-HSD1, PDB ID: 1XU7, UniProt ID: P28845)

For each target, the co-crystal ligand, along with high-affinity and low-affinity GtoPdb ligands, were included as controls.

### Receptor and Ligand Preparation
Protein structures were retrieved from the RCSB Protein Data Bank (PDB). Receptor preparation for docking was performed using Meeko (version 0.7.1) with PDBFixer (version 1.12.0) for initial repair steps. Missing atoms and hydrogens were added, and protonation states were assigned at pH 7.4. Unresolved residue segments were recorded but not modeled. Heterogens were excluded from PDBFixer and restored unchanged after protein repair. A restrained local OpenMM minimization was performed for relaxation, with parameters and energies recorded in the repair audit.

A specific deviation occurred during the preparation of the ERbeta receptor (PDB ID: 1L2J). PDBFixer output for 1L2J created a residue valence conflict in Meeko 0.7.1. To proceed, observed chain-A coordinates were passed directly to Meeko, and Meeko was allowed to exclude observed residues that could not match a complete chemical template. Consequently, missing atoms and unresolved segments were not modeled for ERbeta.

Ligand 3D structures were generated using RDKit ETKDGv3 and prepared for AutoDock Vina using Meeko.

### Molecular Docking
Molecular docking was performed using AutoDock Vina (version 1.2.7) with the Vina scoring function. Docking boxes were derived from the co-crystallized ligand in each PDB structure, with a padding of 5 Å. An exhaustiveness of 4 and 3 poses were requested per ligand. A single seed (20260717) was used for all docking calculations. Each receptor batch of 26 ligands was docked serially on one CPU.

### Data Analysis
For each target, the Vina docking scores (kcal/mol) were obtained. Ligand efficiency (LE) was calculated as the absolute Vina score magnitude divided by the number of heavy atoms. Rank concordance between the Vina and reference MOE docking results was quantified using the Spearman rank correlation coefficient. Additionally, the overlap of the top-5 ranked ligands from both methods was assessed using the Jaccard index. The Vina and MOE scores were not compared for magnitude, only for their relative rankings.

### Execution Provenance and Cost
The computational experiments were executed on Google Cloud Run Jobs, identified as a server tool fallback runtime, as the preferred Google Colab runtime was not used for the docking calculations. The container digest for the `terpedia-cheminformatics` software was `sha256:95e3067ff1592f17200b72852e9f5e92ef9f1a782f1a960259f5c9b096c1eec6`. The total estimated list price for all docking runs was 0.0387 USD.

## Results

A total of 130 Vina docking observations were successfully generated across the five targets. These observations included 23 terpenes and 3 control ligands for each of the five targets. Comparisons were made for 129 shared ligand-target pairs between the newly computed Vina docking results and the reference MOE docking results.

The per-target Spearman rank correlation coefficients between Vina and MOE docking results were:
*   **11beta-HSD1:** 0.559 (moderate)
*   **AR:** 0.2752 (low)
*   **CB2:** 0.292 (low)
*   **CYP19A1:** 0.4435 (modest)
*   **ERbeta:** 0.4161 (modest)

The top-5 ligand overlaps, quantified by the Jaccard index, were:
*   **11beta-HSD1:** 0.6667 (4 out of 6 ligands in common)
    *   Current Top 5: enoxolone, CHAPS, AZD4017, alpha-bisabolol, bornyl acetate
    *   Reference Top 5: CHAPS, AZD4017, enoxolone, nerolidol, alpha-bisabolol
*   **AR:** 0.25 (2 out of 8 ligands in common)
    *   Current Top 5: flutamide, alpha-bisabolol, nerolidol, R-carvone, perillyl alcohol
    *   Reference Top 5: CHEMBL2177238, flutamide, testosterone propionate, nerolidol, bornyl acetate
*   **CB2:** 0.6667 (4 out of 6 ligands in common)
    *   Current Top 5: E3R, AM10257, ibipinabant, alpha-bisabolol, bornyl acetate
    *   Reference Top 5: AM10257, ibipinabant, E3R, nerolidol, alpha-bisabolol
*   **CYP19A1:** 0.6667 (4 out of 6 ligands in common)
    *   Current Top 5: androstenedione, testolactone, letrozole, alpha-bisabolol, nerolidol
    *   Reference Top 5: bornyl acetate, nerolidol, alpha-bisabolol, letrozole, androstenedione
*   **ERbeta:** 1.0 (5 out of 5 ligands in common)
    *   Current Top 5: R,R-THC, 4-hydroxytamoxifen, ospemifene, alpha-bisabolol, nerolidol
    *   Reference Top 5: 4-hydroxytamoxifen, ospemifene, R,R-THC, nerolidol, alpha-bisabolol

The ligand efficiency for each docked ligand was calculated as the absolute Vina score magnitude divided by its heavy-atom count. For example, for AZD4017 at 11beta-HSD1, with a Vina score of -8.956 kcal/mol and 29 heavy atoms, the ligand efficiency was 0.3088. For CHAPS at 11beta-HSD1, with a Vina score of -9.108 kcal/mol and 42 heavy atoms, the ligand efficiency was 0.2169.

## Discussion

This study aimed to provide a reproducible sensitivity analysis of a prior molecular docking investigation of terpenes by replacing a proprietary docking engine (MOE) with an open-source alternative (AutoDock Vina). The results indicate varying degrees of rank concordance between the two docking methodologies across the five selected protein targets.

The Spearman rank correlation coefficients ranged from low (AR: 0.2752, CB2: 0.292) to modest (ERbeta: 0.4161, CYP19A1: 0.4435) and moderate (11beta-HSD1: 0.559). These correlations suggest that while there is some agreement in ligand ranking between MOE and Vina, the extent of this agreement is highly target-dependent. It is important to note that these correlations reflect the consistency of ranking, not the numerical equivalence of the docking scores, which are derived from different scoring functions and are not directly comparable in magnitude.

The top-5 ligand overlap analysis provides an exploratory measure of cross-engine agreement for the highest-ranked ligands. While ERbeta showed a perfect Jaccard index of 1.0, indicating identical top-5 ligands, this finding must be interpreted with caution due to a known deviation in the ERbeta receptor preparation protocol. For other targets, the overlaps were also substantial (e.g., 11beta-HSD1, CB2, CYP19A1 all at 0.6667), suggesting that despite differences in scoring functions and preparation methods, the two engines can identify some common high-ranking ligands. However, top-five overlap is an exploratory cross-engine agreement and does not establish binding, affinity, predictive accuracy, or experimental validity. The discrepancies observed, particularly the lower overlap for AR, highlight the sensitivity of docking results to the specific algorithms and parameters employed.

The calculation of ligand efficiency provides a normalized metric for comparing ligand potency irrespective of size. This metric, defined as the absolute Vina score magnitude divided by the heavy-atom count, allows for a more equitable comparison of ligands, especially when dealing with a diverse set like terpenes.

Overall, this reproduction study demonstrates that an open-source workflow can provide valuable insights into the robustness and sensitivity of computational docking predictions. The observed rank concordances, ranging from low to moderate, underscore the importance of considering methodological choices in molecular docking studies and the need for careful interpretation of results, especially when comparing across different software platforms.

## Limitations

Several limitations are present in this reproduction study. Firstly, the unavailability of row-level outputs from the four target-fishing engines used in the original study prevented the recomputation of target-fishing consensus, TF+/TF- enrichment, and tool ablation analyses. This restricted the scope of the reproduction to only the molecular docking component.

Secondly, a specific deviation occurred during the preparation of the ERbeta receptor (PDB ID: 1L2J). Due to a residue valence conflict with PDBFixer output in Meeko 0.7.1, observed chain-A coordinates were passed directly to Meeko, and Meeko was allowed to exclude incomplete residues that could not match a chemical template. Consequently, missing atoms and unresolved segments were not modeled for ERbeta, which may influence the docking results for this specific target.

Finally, while this study provides a sensitivity analysis of ligand ranking between two different docking engines, the results do not establish experimental binding, affinity, or predictive accuracy. Molecular docking is a computational prediction method, and its outputs require experimental validation.

## Conclusions

This study successfully reproduced the molecular docking component of a previous terpene investigation using an open-source RDKit/Meeko/AutoDock Vina workflow. A total of 130 Vina docking observations were generated, and 129 shared ligand-target pairs were compared with reference MOE results. Spearman rank correlations between Vina and MOE docking results ranged from low (AR: 0.2752, CB2: 0.292) to modest (ERbeta: 0.4161, CYP19A1: 0.4435) and moderate (11beta-HSD1: 0.559). Top-5 ligand overlaps also varied, with ERbeta showing a perfect overlap (1.0, 5/5 ligands) despite a receptor preparation deviation. This work demonstrates the utility of open-source tools for performing reproducible sensitivity analyses in computational drug discovery, highlighting the influence of methodological choices on predicted binding profiles.

## Data Availability

The research data, including the executed low-compute study records and the Google Colab notebook (revision 4), are available in a private GitHub repository and a content-addressed Google Colab notebook.

*   **GitHub Repository:** `https://github.com/Terpedia/terpedia-reproduce-the-research-in-the-attached-reference-manusc-be6fa4a4`
*   **Colab Notebook URL:** `https://colab.research.google.com/github/Terpedia/terpedia-reproduce-the-research-in-the-attached-reference-manusc-be6fa4a4/blob/main/notebooks/study_be6fa4a4-c9eb-45e6-916e-aec999992e4f-reproducibility-r4.ipynb`
*   **Research Snapshot SHA256:** `4c08c84feac979e8ba665156a6b8f42c10c7ef741f3912a89b3b5b0ea2e2ee33`
*   **Notebook Content SHA256:** `5d5787b35deade1dee2f6d8d7a895947458e15be361de50bc39890e53ffbfee0`
*   **Study Bundle SHA256:** `95c61eae810bc51d6721dfc29ea7703807608bedeb43a202a02156271f580f79`
*   **Total Estimated Compute Cost:** 0.0387 USD.
*   **Execution Platform:** Google Cloud Run Jobs (server tool fallback runtime) for docking calculations.

## AI Disclosure

This manuscript was authored and revised by the Terpedia Editor Agent, GPT-420, an AI assistant. The agent utilized supplied study records, a research notebook, and specific editorial instructions to generate and refine the content, ensuring adherence to scientific reporting standards and preservation of data provenance.

## References
[1] Placeholder citation 1
[2] Placeholder citation 2
[3] Placeholder citation 3