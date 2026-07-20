# Evidence-qualified reconstruction of the *Cannabis sativa* metabolic network from Terpedia

## Abstract

We reconstructed a provenance-linked draft metabolic network for *Cannabis
sativa* in Terpedia and evaluated whether its reaction reactants are reachable
from an explicit plant defined-medium proxy. The reconstruction uses the UniProt
reference proteome UP000583929 (30,304 proteins; Jamaican Lion 4 male) joined to
Rhea reactions and ChEBI chemical participants. The resulting catalog contains
1,086 Rhea master reactions and 1,245 catalog ChEBI participants, of which 831
unique metabolites occur as reaction reactants in the reachability graph. Only
16 reactants (1.93%) are reachable under curated physiological directions, 24
(2.89%) under a hybrid direction mode, and 26 (3.13%) under a fully reversible
structural upper bound. The analysis produces 70 unresolved exact-EC reaction
hypotheses, 369 missing-producer hypotheses, and 1,319 blocked-known-reaction
hypotheses; 1,158 blocked hypotheses retain candidate enzyme/protein support.
The catalog has 3,564 proteins with exact EC annotations, 695 exact EC numbers
mapped to Rhea, 1,032 reactions with directional reaction SMILES, 164 reactions
with physiological direction evidence, and 1,487 transporter candidates. This
is not a complete Cannabis flux model: compartments, exchanges, biomass,
thermodynamic constraints, and complete gene–protein–reaction logic remain
uncurated. Terpedia therefore provides a reproducible, evidence-qualified
research substrate and a prioritized hypothesis agenda, not proof of every
Cannabis pathway or cannabinoid biosynthetic flux.

**Keywords:** *Cannabis sativa*; metabolic reconstruction; pathway inference;
Rhea; ChEBI; enzyme annotation; cannabinoid biology; Terpedia

## 1. Introduction

Cannabis research spans a chemically rich metabolome, cultivar variation, and a
genome whose enzyme annotations are distributed across heterogeneous databases.
The presence of a compound, a pathway association, or an enzyme annotation is
not equivalent to a complete organismal metabolic model. A useful reconstruction
must retain those distinctions and expose where the graph cannot bootstrap from
defined-medium inputs.

This study applies the Terpedia pathway-inference workflow previously used for
fungal reconstruction to *C. sativa*. We ask three questions: (i) how many
reference-proteome enzyme annotations resolve to standardized reactions; (ii)
what fraction of unique reaction reactants is reachable from a plant defined-
medium proxy; and (iii) which reaction, enzyme, and gene hypotheses should be
prioritized to close the network.

## 2. Materials and methods

### 2.1 Reference proteome and reaction projection

The reference proteome was UniProt UP000583929, *C. sativa* cv. Jamaican Lion 4
(male), containing 30,304 protein entries and associated with the JAATIQ010000000
assembly. Exact and partial EC annotations, direct catalytic activities,
physiological direction evidence, candidate compartments, and transport-related
annotations were parsed from the UniProt tabular export. Exact EC numbers were
joined to Rhea master reactions. Rhea equations were parsed into ChEBI reactant
and product entities; formulas, structures, elemental balance, and formal charge
were retained where available.

This projection is an annotation catalog, not a manually validated Cannabis
model. Exact-EC-to-Rhea joins are explicitly distinguished from direct
biochemical assays, homology, expression, and metabolomics observation.

### 2.2 Hypergraph reachability

Each reaction is a directed hyperedge. An operation fires only when every
metabolite on its input side is reachable. Products are added until a fixed
point. The denominator is every unique ChEBI entity appearing as a reactant in
the reaction graph, yielding 831 reactants; the larger catalog participant count
(1,245) includes ChEBI entities retained outside the executable reaction sides.

The seed set is a glucose-based plant defined-medium proxy containing glucose,
ammonium, nitrate, phosphate, sulfate, potassium, sodium, magnesium, calcium,
iron, molybdate, water, proton, oxygen, and carbon dioxide. Eighteen network
entities resolve to these seeds. Zinc, manganese, borate, and biotin are retained
as requested ingredients but are absent from the reaction-participant catalog.
The proxy is computational and is not asserted to be a tissue-culture recipe for
Cannabis.

Three direction modes were evaluated: curated physiological directions only;
curated directions plus unknown-direction reactions treated as reversible; and
all reactions reversible as a structural upper bound. Stoichiometric
coefficients were retained but not interpreted as flux quantities. No exchange,
biomass, or compartment constraints were invented.

### 2.3 Hypothesis generation

Three hypothesis classes were generated. Exact-EC proteins without a projected
`catalyzes` edge were grouped by EC number as unresolved reaction hypotheses.
Unreachable metabolites with no producer reaction in the included Rhea catalog
were recorded as missing-producer hypotheses. Existing reactions consuming
unreachable metabolites were recorded as blocked-known-reaction hypotheses,
retaining the equation and available protein, locus, gene, and EC support.
Missing reaction chemistry was never fabricated.

### 2.4 Reproducibility

The analysis was generated in the Terpedia Knowledge repository and published
with the companion Colab notebook in this research repository. The notebook
clones Terpedia Knowledge, loads the versioned Cannabis network and reports,
reconstructs the reaction sides, reruns the fixed-point closure, and asserts
agreement with the stored results.

```bash
npm run build:cannabis-sativa-network
npm run build:cannabis-sativa-reachability
npm run build:cannabis-sativa-hypotheses
npm test
```

## 3. Results

### 3.1 Annotation and catalog scale

The reference proteome contains 3,564 proteins with exact EC annotations
(11.76%), representing 765 unique exact EC numbers and 41 partial EC numbers.
Six hundred ninety-five exact EC numbers map to Rhea (90.85%), resolving 2,885
proteins in the source annotation calculation. The projected catalog contains
1,086 Rhea master reactions, 1,032 reactions with directional reaction SMILES,
and 1,245 catalog ChEBI participants. Of these reactions, 676 are element-
balanced and 801 are charge-balanced among auditable equations; 410 and 285,
respectively, are not auditable because of generalized or incomplete chemical
representations.

The reaction graph has 21 bipartite connected components. Its largest component
contains 1,063 reactions and 1,195 metabolites, while 634 metabolites are
boundary participants. These are structural catalog statistics and do not imply
cellular flux.

### 3.2 Defined-medium reachability

| Direction mode | Reachable reactants | Reactant denominator | Reachability |
|---|---:|---:|---:|
| Curated directions only | 16 | 831 | 1.93% |
| Curated directions + unknown reversible | 24 | 831 | 2.89% |
| Fully reversible upper bound | 26 | 831 | 3.13% |

Relaxing directionality adds only two percentage points of reactant closure in
absolute terms. The low upper bound indicates that missing producer chemistry,
unreachable precursor chains, and unseeded intracellular currencies dominate the
current gap. ATP, NAD(P)H, CoA, and amino acids were not silently seeded because
doing so would change the defined-medium question.

### 3.3 Reaction, enzyme, and gene hypotheses

The hypothesis layer identifies 684 exact-EC proteins without a projected
`catalyzes` edge, grouped into 70 unresolved exact-EC hypotheses. A further 478
proteins have partial EC annotations without an included Rhea mapping. There are
369 metabolites with no producer reaction in the catalog. Existing reactions
produce 1,319 blocked-known-reaction hypotheses; 1,158 retain at least one
candidate enzyme/protein association.

These hypotheses are curation targets, not newly asserted Cannabis biochemistry.
The highest-value next steps are authoritative EC-to-Rhea resolution, producer
reaction discovery for central precursor gaps, and validation of candidate
enzyme–locus links with reciprocal orthology, tissue expression, and biochemical
assays.

### 3.4 Compartments, transport, and directionality

The catalog contains 1,487 transporter candidates, 1,210 with transmembrane
topology support and 801 with annotation-specificity support. Twenty-four Rhea
reactions have explicit in/out notation, but none are promoted to curated
Cannabis transport or exchange reactions. Candidate compartment evidence is
available for 728 reactions (67.03%), while only 164 reactions (15.10%) have
physiological direction evidence and 11 have experimental direction evidence.

### 3.5 Cannabis metabolome context

Terpedia’s CannabisDB layer provides compound identity, cultivar concentration,
protein-association, and pathway-association records. Those relations are
valuable metabolome and target context, but they do not supply stoichiometric
reaction equations or establish biosynthetic flux. Accordingly, this network
analysis does not report a cannabinoid pathway “coverage percentage.” A future
version should define an explicit, source-grounded cannabinoid reaction panel
before calculating such a metric.

## 4. Critical review and threats to validity

The headline percentage is conditional on the denominator, seed set, and graph
semantics. Using all participants instead of reactants, or seeding ATP and
NAD(P)H, would produce different values. The plant medium proxy is not a
validated Cannabis tissue-culture formulation. Direction evidence is annotation
support rather than thermodynamic feasibility. Exact-EC, projected-Rhea,
candidate-enzyme, transporter, and metabolomics counts use different entities
and denominators and must not be combined into a global score.

The 684-protein hypothesis count is the number of exact-EC protein entities with
no projected `catalyzes` edge in this Terpedia bundle. It is not the arithmetic
complement of the separate 2,885-protein source annotation-resolution value.
This distinction is necessary because the joins answer different questions.

The Colab notebook independently reruns the graph closure and checks stored
counts, but it does not independently validate upstream UniProt, Rhea, ChEBI,
or CannabisDB records. A peer-review release should pin a commit and archive
all upstream source files and release metadata.

## 5. Discussion

The Cannabis catalog is larger than the current *P. cubensis* catalog in both
proteome and reaction counts, but its reachability remains low. This is not
evidence that Cannabis lacks metabolism; it shows that a reaction catalog built
from available enzyme annotations does not yet close from a defined medium.
The reversible upper bound of 3.13% is a useful diagnostic of missing precursor
and cofactor closure rather than a biological activity estimate.

The practical reconstruction agenda is to resolve the 70 exact-EC hypotheses,
add experimentally supported producer reactions for the 369 no-producer
metabolites, and then add explicit compartments, transport, exchanges, biomass,
and Boolean gene–protein–reaction rules. CannabisDB compound and cultivar data
can then be used as orthogonal metabolome context, while curated cannabinoid
reaction panels can be evaluated without conflating pathway association with
stoichiometric reaction coverage.

## 6. Limitations

The reference proteome represents one male cultivar and is not population-
complete. Automated protein annotations may be incomplete or overly broad.
Generalized Rhea participants limit elemental and charge auditing. Candidate
compartments and transport annotations are not model assignments. CannabisDB
pathway and target relations do not establish reaction direction or flux. No
biomass objective, exchange system, thermodynamic model, or inferential
differential-expression layer is included.

## 7. Conclusions

Terpedia now supports a reproducible *C. sativa* draft metabolic reconstruction
with explicit reachability and reaction/enzyme/gene hypotheses. Only 3.13% of
unique reaction reactants are reachable under the fully reversible upper bound
from the plant defined-medium proxy. The network should be treated as an
auditable research substrate and prioritized curation agenda, not as a complete
Cannabis metabolic or flux-balance model.

## Data and code availability

The computational source and generated artifacts are in the Terpedia Knowledge
repository under the Cannabis network, reachability, and hypothesis scripts and
reports. The companion notebook is published in this repository at
`notebooks/cannabis_sativa_metabolism_colab.ipynb`.

## References

1. UniProt *Cannabis sativa* reference proteome UP000583929:
   <https://www.uniprot.org/proteomes/UP000583929>
2. Rhea reaction database: <https://www.rhea-db.org/>
3. ChEBI chemical ontology: <https://www.ebi.ac.uk/chebi/>
4. Cannabis Database (CannabisDB) source relations:
   <https://cannabisdatabase.ca/>
5. Terpedia Knowledge source repository:
   <https://github.com/Terpedia/terpedia-knowledge>
