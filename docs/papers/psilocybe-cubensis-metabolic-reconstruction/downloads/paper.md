# Evidence-qualified reconstruction of the *Psilocybe cubensis* metabolic network from Terpedia

## Abstract

We reconstructed a provenance-linked draft metabolic network for *Psilocybe
cubensis* by joining the UniProt UP000664032 reference proteome to Rhea
reactions and ChEBI chemical participants. The catalog contains 826 biochemical
reactions and 964 unique ChEBI participants. We defined completeness as the
fraction of unique reaction reactants reachable from a glucose-based fungal
defined-medium proxy using a forward reaction-hypergraph closure. Twelve of 617
reactants (1.94%) are reachable under curated physiological directions; 18
(2.92%) under a hybrid mode; and 18 (2.92%) under the fully reversible
structural upper bound. The reconstruction generates 63 unresolved exact-EC
reaction hypotheses, 288 missing-producer hypotheses, and 1,007 blocked-known-
reaction hypotheses, with 855 retaining candidate enzyme/protein support. All
eight expected Rhea reactions in the curated psilocybin subnetwork are
represented, but the global network is not a complete flux model. Terpedia
therefore provides an auditable reaction catalog and prioritized hypothesis
agenda rather than proof of complete *in vivo* metabolism.

**Keywords:** *Psilocybe cubensis*; metabolic reconstruction; Rhea; ChEBI;
psilocybin; pathway inference; Terpedia

## 1. Introduction

Genome annotation, reaction representation, metabolite detection, and pathway
activity are different claims. This study uses Terpedia to preserve those
distinctions while measuring whether the reaction catalog can bootstrap from a
defined-medium seed set. The aims were to quantify reaction resolution, measure
reactant reachability, and generate evidence-qualified reaction, enzyme, and
gene hypotheses for closing the network.

## 2. Methods

The UniProt UP000664032 reference proteome was joined to Rhea through exact EC
annotations. Rhea equations were parsed into explicit `has_reactant` and
`has_product` ChEBI edges, with molecular formula, structure, balance, candidate
compartment, transporter, and physiological-direction evidence retained.

The denominator was every unique ChEBI entity appearing as a reaction reactant
(617). The seed set was a glucose-based fungal defined-medium proxy containing
glucose, ammonium, nitrate, phosphate, sulfate, chloride, calcium, iron, copper,
biotin, water, proton, oxygen, and carbon dioxide. Three modes were evaluated:
curated directions only; curated directions plus unknown reactions reversible;
and all reactions reversible. Stoichiometry was retained but not interpreted as
flux, and no biomass, exchange, thermodynamic, or compartment constraints were
invented.

Exact-EC proteins without a projected `catalyzes` edge were grouped as
unresolved reaction hypotheses. Unreachable metabolites with no producer
reaction were recorded as missing-producer hypotheses. Reactions consuming
unreachable metabolites were recorded as blocked-known-reaction hypotheses,
retaining existing equations and enzyme/locus/gene support. No missing reaction
equations were fabricated.

## 3. Results

### 3.1 Catalog

The proteome contains 13,113 proteins, 1,107 with exact EC annotations (8.44%).
The reaction catalog contains 826 Rhea master reactions, 772 reactions with
directional reaction SMILES, and 964 unique ChEBI participants. All equations
parsed; 489 were element-auditable and 582 charge-auditable. The participant
graph has 20 connected components; its largest contains 806 reactions and 937
metabolites.

### 3.2 Reachability

| Direction mode | Reachable reactants | Denominator | Reachability |
|---|---:|---:|---:|
| Curated directions only | 12 | 617 | 1.94% |
| Curated + unknown reversible | 18 | 617 | 2.92% |
| Fully reversible upper bound | 18 | 617 | 2.92% |

The identical hybrid and reversible values indicate that missing entry chemistry
and precursor/cofactor closure dominate the gap rather than directionality alone.

### 3.3 Hypotheses and psilocybin subnetwork

The hypothesis layer contains 63 unresolved exact-EC hypotheses, 288 metabolites
with no producer reaction, and 1,007 blocked-known-reaction hypotheses. Eight of
eight expected Rhea reactions in the curated psilocybin subnetwork are
represented. This is reaction representation coverage, not evidence of pathway
flux in every tissue or strain.

### 3.4 Model readiness

Two public RNA-seq runs contribute descriptive reaction-expression links. There
are 447 transporter candidates, but zero curated transport or exchange reactions.
Only 144 reactions have physiological direction evidence. The network lacks a
biomass objective, complete Boolean gene–protein–reaction rules, curated
compartments, thermodynamic constraints, and a validated exchange system.

## 4. Critical review and limitations

The reachability percentage is conditional on its denominator, seed formulation,
and direction assumptions. It is not a global completeness score and should not
be compared directly with genome BUSCO completeness or metabolomics feature
counts. The fungal defined-medium formulation is a computational proxy, not a
*P. cubensis*-validated growth medium. Exact-EC joins and RNA expression do not
establish enzyme activity; metabolite detection does not establish biosynthetic
origin. Generalized Rhea participants limit elemental and charge auditing.

The correct next steps are to resolve exact EC-to-Rhea gaps, add authoritative
producer reactions, run explicitly labeled cofactor sensitivity scenarios, and
curate compartments, transport, exchange, biomass, and gene–protein–reaction
logic.

## 5. Conclusion

Terpedia provides a reproducible, evidence-qualified *P. cubensis* metabolic
reconstruction. Only 2.92% of unique reaction reactants are reachable under the
reversible upper bound. The artifact is best used as a transparent curation and
experimental agenda, not as a finished organismal or flux-balance model.

## Data and code availability

Source artifacts and reconstruction scripts are in the [Terpedia Knowledge
repository](https://github.com/Terpedia/terpedia-knowledge). The companion
Colab notebook is published with this paper.

## References

1. UniProt reference proteome UP000664032: <https://www.uniprot.org/proteomes/UP000664032>
2. Rhea reaction database: <https://www.rhea-db.org/>
3. ChEBI chemical ontology: <https://www.ebi.ac.uk/chebi/>
4. A draft reference assembly of the *Psilocybe cubensis* genome:
   <https://doi.org/10.12688/f1000research.51613.2>
5. Enzymatic Synthesis of Psilocybin: <https://doi.org/10.1002/anie.201705489>
