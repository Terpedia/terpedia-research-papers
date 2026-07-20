from typing import Literal

from pydantic import BaseModel, Field, model_validator


StructureFormat = Literal["auto", "smiles", "inchi", "molblock", "name"]
DockingEngine = Literal["gnina", "vina"]
CnnScoringMode = Literal["rescore", "refinement"]
CnnModel = Literal["default", "fast"]


class StructureInput(BaseModel):
    structure: str = Field(min_length=1, max_length=100_000)
    format: StructureFormat = "auto"
    label: str | None = Field(default=None, max_length=200)


class AnalyzeRequest(StructureInput):
    depict: bool = True
    conformers: int = Field(default=0, ge=0, le=20)
    include_alerts: bool = True


class CompareRequest(BaseModel):
    molecules: list[StructureInput] = Field(min_length=2, max_length=25)
    fingerprint: Literal["morgan", "rdkit", "maccs"] = "morgan"
    radius: int = Field(default=2, ge=1, le=4)
    n_bits: int = Field(default=2048, ge=256, le=8192)
    include_mcs: bool = True
    depict: bool = True


class SubstructureRequest(BaseModel):
    query: StructureInput
    targets: list[StructureInput] = Field(min_length=1, max_length=100)
    use_chirality: bool = True
    depict: bool = True


class ReactionRequest(BaseModel):
    reaction: str = Field(min_length=3, max_length=20_000)
    format: Literal["smarts", "smiles"] = "smiles"
    depict: bool = True


class ResolveRequest(BaseModel):
    name: str = Field(min_length=1, max_length=300)
    analyze: bool = True
    depict: bool = True


class PharmacophoreRequest(StructureInput):
    conformers: int = Field(default=1, ge=1, le=10)


class MetabolismRequest(StructureInput):
    max_depth: int = Field(default=2, ge=1, le=3)
    max_products: int = Field(default=24, ge=1, le=60)
    phases: list[Literal["I", "II"]] = Field(default_factory=lambda: ["I", "II"], min_length=1, max_length=2)
    depict: bool = True


class Vector3(BaseModel):
    x: float = Field(ge=-1000, le=1000)
    y: float = Field(ge=-1000, le=1000)
    z: float = Field(ge=-1000, le=1000)

    def values(self) -> list[float]:
        return [self.x, self.y, self.z]


class DockingBox(BaseModel):
    center: Vector3
    size: Vector3

    @model_validator(mode="after")
    def validate_size(self):
        dimensions = self.size.values()
        if any(value < 8 or value > 40 for value in dimensions):
            raise ValueError("docking box dimensions must be between 8 and 40 Angstrom")
        if dimensions[0] * dimensions[1] * dimensions[2] > 27_000:
            raise ValueError("docking box volume cannot exceed 27,000 cubic Angstrom")
        return self


class ReferenceLigand(BaseModel):
    residue_name: str = Field(min_length=1, max_length=8, pattern=r"^[A-Za-z0-9]+$")
    chain: str | None = Field(default=None, max_length=4)
    residue_number: int | None = Field(default=None, ge=-9999, le=99999)
    padding: float = Field(default=5.0, ge=2.0, le=12.0)


class ReceptorInput(BaseModel):
    pdb_id: str | None = Field(default=None, pattern=r"^[0-9][A-Za-z0-9]{3}$")
    pdb_text: str | None = Field(default=None, min_length=20, max_length=5_000_000)
    pdbqt_text: str | None = Field(default=None, min_length=20, max_length=5_000_000)
    chains: list[str] = Field(default_factory=list, max_length=8)
    keep_heterogens: list[str] = Field(default_factory=list, max_length=16)

    @model_validator(mode="after")
    def exactly_one_source(self):
        sources = [self.pdb_id, self.pdb_text, self.pdbqt_text]
        if sum(value is not None for value in sources) != 1:
            raise ValueError("provide exactly one receptor source: pdb_id, pdb_text, or pdbqt_text")
        self.pdb_id = self.pdb_id.upper() if self.pdb_id else None
        self.chains = list(dict.fromkeys(chain.strip() for chain in self.chains if chain.strip()))
        self.keep_heterogens = list(dict.fromkeys(name.upper() for name in self.keep_heterogens))
        return self


class DockingRequest(BaseModel):
    ligand: StructureInput
    receptor: ReceptorInput
    box: DockingBox | None = None
    reference_ligand: ReferenceLigand | None = None
    engine: DockingEngine = "gnina"
    scoring: Literal["vina", "vinardo"] = "vina"
    cnn_scoring: CnnScoringMode = "rescore"
    cnn_model: CnnModel = "default"
    exhaustiveness: int = Field(default=8, ge=1, le=32)
    poses: int = Field(default=5, ge=1, le=10)
    energy_range: float = Field(default=4.0, gt=0, le=10)
    seed: int = Field(default=20260717, ge=1, le=2_147_483_647)

    @model_validator(mode="after")
    def require_site(self):
        if self.box is None and self.reference_ligand is None and not self.receptor.pdb_id:
            raise ValueError("provide a docking box or reference ligand for non-PDB receptor input")
        if self.receptor.pdbqt_text and self.box is None:
            raise ValueError("prepared PDBQT receptors require an explicit docking box")
        return self


class BatchDockingRequest(BaseModel):
    ligands: list[StructureInput] = Field(min_length=1, max_length=30)
    receptor: ReceptorInput
    box: DockingBox | None = None
    reference_ligand: ReferenceLigand | None = None
    engine: DockingEngine = "gnina"
    scoring: Literal["vina", "vinardo"] = "vina"
    cnn_scoring: CnnScoringMode = "rescore"
    cnn_model: CnnModel = "fast"
    exhaustiveness: int = Field(default=4, ge=1, le=16)
    poses: int = Field(default=3, ge=1, le=5)
    energy_range: float = Field(default=4.0, gt=0, le=10)
    seed: int = Field(default=20260717, ge=1, le=2_147_483_647)
    include_pose_artifacts: bool = False

    @model_validator(mode="after")
    def require_site(self):
        if self.box is None and self.reference_ligand is None and not self.receptor.pdb_id:
            raise ValueError("provide a docking box or reference ligand for non-PDB receptor input")
        if self.receptor.pdbqt_text and self.box is None:
            raise ValueError("prepared PDBQT receptors require an explicit docking box")
        return self


class HqBatchDockingRequest(BatchDockingRequest):
    engine: DockingEngine = "gnina"
    cnn_scoring: CnnScoringMode = "refinement"
    cnn_model: CnnModel = "default"
    exhaustiveness: int = Field(default=32, ge=8, le=64)
    poses: int = Field(default=10, ge=3, le=20)
    seeds: list[int] = Field(
        default_factory=lambda: [20260717, 20260718, 20260719, 20260720, 20260721],
        min_length=3, max_length=7,
    )
    scoring_functions: list[Literal["vina", "vinardo"]] = Field(
        default_factory=lambda: ["vina", "vinardo"], min_length=2, max_length=2
    )
    validation_exhaustiveness: int = Field(default=32, ge=8, le=64)
    redocking_rmsd_threshold: float = Field(default=2.0, gt=0, le=5)
    require_redocking_pass: bool = True

    @model_validator(mode="after")
    def unique_hq_parameters(self):
        self.seeds = list(dict.fromkeys(self.seeds))
        self.scoring_functions = list(dict.fromkeys(self.scoring_functions))
        if len(self.seeds) < 3:
            raise ValueError("HQ docking requires at least three unique seeds")
        if self.engine == "vina" and set(self.scoring_functions) != {"vina", "vinardo"}:
            raise ValueError("HQ docking requires both vina and vinardo scoring")
        return self


class BatchAnalyzeRequest(BaseModel):
    molecules: list[AnalyzeRequest] = Field(min_length=1, max_length=100)

    @model_validator(mode="after")
    def bound_work(self):
        if sum(item.conformers for item in self.molecules) > 100:
            raise ValueError("batch conformer total cannot exceed 100")
        return self


class TargetPredictionRecord(BaseModel):
    ligand: str = Field(min_length=1, max_length=200)
    target: str = Field(min_length=1, max_length=200)
    tool: Literal["SwissTargetPrediction", "SEA", "STarFish", "PPB2"]


class DockingObservation(BaseModel):
    ligand: str = Field(min_length=1, max_length=200)
    target: str = Field(min_length=1, max_length=200)
    score_kcal_mol: float
    heavy_atoms: int = Field(ge=1, le=500)
    ligand_class: Literal["monoterpene", "sesquiterpene", "known_ligand"]
    experimental_support: bool = False


class ExperimentalDockingPrediction(BaseModel):
    ligand: str = Field(min_length=1, max_length=200)
    pubchem_cid: int | None = Field(default=None, ge=1)
    target: str = Field(min_length=1, max_length=200)
    score: float


class ExperimentalValidationRequest(BaseModel):
    predictions: list[ExperimentalDockingPrediction] = Field(min_length=2, max_length=20_000)
    evidence_scope: Literal["direct_binding", "functional"] = "direct_binding"
    lower_is_better: bool = True
    active_threshold_nm: float = Field(default=10_000, gt=0, le=1_000_000)
    score_name: str = Field(default="docking_score", min_length=1, max_length=80)


class StudyReplicationRequest(BaseModel):
    question: str = Field(min_length=1, max_length=4000)
    predictions: list[TargetPredictionRecord] = Field(default_factory=list, max_length=20_000)
    docking: list[DockingObservation] = Field(default_factory=list, max_length=20_000)
    reference_docking: list[DockingObservation] = Field(default_factory=list, max_length=20_000)
    reference_method: str = Field(default="Reference paper docking method", min_length=1, max_length=300)
    consensus_min_tools: int = Field(default=2, ge=1, le=4)
    enrichment_fraction: float = Field(default=0.2, gt=0, le=1)
    include_manifest_rows: bool = False
    compute_profile: Literal["low_compute", "balanced", "high_quality"] = "low_compute"
