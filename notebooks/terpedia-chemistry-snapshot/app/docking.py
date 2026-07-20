import datetime as dt
import functools
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from statistics import mean, median, stdev

import httpx
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolAlign

from .chemistry import ChemistryError, ChemistryUpstreamError, SERVICE_VERSION, identifiers
from .receptor_repair import repair_receptor_pdb


RCSB_FILES = "https://files.rcsb.org/download"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core/entry"
MAX_RECEPTOR_ATOMS = int(os.getenv("DOCKING_MAX_RECEPTOR_ATOMS", "25000"))
DOCKING_CPUS = int(os.getenv("DOCKING_CPUS", "1"))
RECEPTOR_REPAIR_PH = float(os.getenv("DOCKING_RECEPTOR_PH", "7.4"))
GNINA_BIN = os.getenv("GNINA_BIN", "gnina")
GNINA_TIMEOUT_SECONDS = int(os.getenv("GNINA_TIMEOUT_SECONDS", "900"))
GNINA_USE_GPU = os.getenv("GNINA_USE_GPU", "false").lower() == "true"

# Reproducible sites from the attached terpene target-fishing manuscript.
TARGET_PRESETS = {
    "6KPF": {
        "name": "Cannabinoid receptor 2", "chains": ["R"],
        "reference": ("E3R", "R", 401), "keep": [],
        "set_templates": {"R:174": "CYX", "R:179": "CYX"},
        "relaxation_iterations": 200,
    },
    "3V49": {"name": "Androgen receptor", "chains": ["A"], "reference": ("PK0", "A", 950), "keep": []},
    "3EQM": {"name": "Cytochrome P450 19A1", "chains": ["A"], "reference": ("ASD", "A", 601), "keep": ["HEM"]},
    "1L2J": {
        "name": "Estrogen receptor beta", "chains": ["A"], "reference": ("ETC", "A", 600), "keep": [],
        "pdbfixer": False,
        "allow_bad_res": True,
        "repair_note": "PDBFixer output for 1L2J creates a residue valence conflict in Meeko 0.7.1; use the observed chain-A coordinates and let Meeko exclude incomplete residues that cannot match a chemical template.",
    },
    "1XU7": {"name": "11-beta-hydroxysteroid dehydrogenase 1", "chains": ["A"], "reference": ("CPS", "A", 523), "keep": ["NDP"]},
}


@dataclass
class PreparedReceptor:
    pdbqt: str
    source: dict
    preparation: dict
    box: dict
    reference_ligand: Chem.Mol | None = None
    reference_label: str | None = None


def _pdb_atom(line: str) -> dict | None:
    if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
        return None
    try:
        return {
            "record": line[:6].strip(),
            "name": line[12:16].strip(),
            "altloc": line[16:17].strip(),
            "residue_name": line[17:20].strip().upper(),
            "chain": line[21:22].strip(),
            "residue_number": int(line[22:26]),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "element": (line[76:78].strip() if len(line) >= 78 else ""),
        }
    except (TypeError, ValueError):
        return None


def box_from_reference_ligand(pdb_text: str, reference) -> dict:
    matches = []
    residue_name = reference.residue_name.upper()
    for line in pdb_text.splitlines():
        atom = _pdb_atom(line)
        if not atom or atom["record"] != "HETATM" or atom["residue_name"] != residue_name:
            continue
        if reference.chain is not None and atom["chain"] != reference.chain:
            continue
        if reference.residue_number is not None and atom["residue_number"] != reference.residue_number:
            continue
        if atom["element"].upper() == "H" or atom["name"].upper().startswith("H"):
            continue
        matches.append(atom)
    if not matches:
        identity = residue_name
        if reference.chain is not None:
            identity += f" chain {reference.chain}"
        if reference.residue_number is not None:
            identity += f" residue {reference.residue_number}"
        raise ChemistryError(f"reference ligand {identity} was not found in the receptor")

    axes = [[atom[key] for atom in matches] for key in ("x", "y", "z")]
    center = [(min(values) + max(values)) / 2 for values in axes]
    dimensions = [max(16.0, max(values) - min(values) + (2 * reference.padding)) for values in axes]
    dimensions = [min(40.0, value) for value in dimensions]
    return {
        "center": {axis: round(value, 3) for axis, value in zip(("x", "y", "z"), center)},
        "size": {axis: round(value, 3) for axis, value in zip(("x", "y", "z"), dimensions)},
        "derived_from": {
            "residue_name": residue_name,
            "chain": reference.chain,
            "residue_number": reference.residue_number,
            "heavy_atoms": len(matches),
            "padding_angstrom": reference.padding,
        },
    }


def extract_reference_ligand(pdb_text: str, reference) -> Chem.Mol:
    selected_serials = set()
    lines = []
    for line in pdb_text.splitlines():
        atom = _pdb_atom(line)
        if not atom or atom["record"] != "HETATM":
            continue
        if atom["residue_name"] != reference.residue_name.upper():
            continue
        if reference.chain is not None and atom["chain"] != reference.chain:
            continue
        if reference.residue_number is not None and atom["residue_number"] != reference.residue_number:
            continue
        lines.append(line)
        try:
            selected_serials.add(int(line[6:11]))
        except ValueError:
            pass
    for line in pdb_text.splitlines():
        if not line.startswith("CONECT"):
            continue
        serials = []
        for offset in range(6, len(line), 5):
            try:
                serials.append(int(line[offset:offset + 5]))
            except ValueError:
                pass
        if serials and serials[0] in selected_serials:
            lines.append(line)
    lines.append("END")
    mol = Chem.MolFromPDBBlock("\n".join(lines) + "\n", sanitize=True, removeHs=True, proximityBonding=True)
    if mol is None or mol.GetNumHeavyAtoms() < 3:
        raise ChemistryError(f"could not reconstruct co-crystal ligand {reference.residue_name} for redocking")
    return mol


def filter_receptor_pdb(pdb_text: str, chains: list[str], keep_heterogens: list[str], reference) -> tuple[str, int]:
    output = []
    atom_count = 0
    chain_set = set(chains)
    keep_set = {name.upper() for name in keep_heterogens}
    for line in pdb_text.splitlines():
        atom = _pdb_atom(line)
        if atom:
            if chain_set and atom["chain"] not in chain_set:
                continue
            if atom["record"] == "HETATM":
                is_reference = (
                    reference is not None
                    and atom["residue_name"] == reference.residue_name.upper()
                    and (reference.chain is None or atom["chain"] == reference.chain)
                    and (reference.residue_number is None or atom["residue_number"] == reference.residue_number)
                )
                if is_reference or atom["residue_name"] not in keep_set:
                    continue
            output.append(line)
            atom_count += 1
        elif line.startswith(("TER", "END")):
            output.append(line)
    if not atom_count:
        raise ChemistryError("receptor filtering removed every atom")
    if atom_count > MAX_RECEPTOR_ATOMS:
        raise ChemistryError(f"receptor has {atom_count} atoms; maximum is {MAX_RECEPTOR_ATOMS}")
    if not output or output[-1] != "END":
        output.append("END")
    return "\n".join(output) + "\n", atom_count


async def fetch_rcsb_receptor(pdb_id: str) -> tuple[str, dict]:
    pdb_id = pdb_id.upper()
    try:
        async with httpx.AsyncClient(timeout=20, follow_redirects=True) as client:
            pdb_response, metadata_response = await __import__("asyncio").gather(
                client.get(f"{RCSB_FILES}/{pdb_id}.pdb"),
                client.get(f"{RCSB_DATA}/{pdb_id}"),
            )
        pdb_response.raise_for_status()
        metadata_response.raise_for_status()
    except (httpx.HTTPError, ValueError) as exc:
        raise ChemistryUpstreamError(f"RCSB could not provide PDB {pdb_id}: {exc}") from exc
    metadata = metadata_response.json()
    title = (metadata.get("struct") or {}).get("title")
    return pdb_response.text, {
        "database": "RCSB Protein Data Bank",
        "pdb_id": pdb_id,
        "title": title,
        "structure_url": f"https://www.rcsb.org/structure/{pdb_id}",
        "coordinate_url": f"{RCSB_FILES}/{pdb_id}.pdb",
        "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }


def prepare_receptor_pdbqt(
    pdb_text: str,
    set_templates: dict[str, str] | None = None,
    allow_bad_res: bool = False,
) -> str:
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer, ResidueChemTemplates

        polymer = Polymer.from_pdb_string(
            pdb_text,
            ResidueChemTemplates.create_from_defaults(),
            MoleculePreparation(),
            set_template=set_templates or {},
            allow_bad_res=allow_bad_res,
            default_altloc="A",
        )
        rigid_pdbqt, flexible = PDBQTWriterLegacy.write_from_polymer(polymer)
        if flexible:
            raise ChemistryError("unexpected flexible residues in rigid receptor preparation")
        if not rigid_pdbqt.strip():
            raise ChemistryError("Meeko produced an empty receptor")
        return rigid_pdbqt
    except ChemistryError:
        raise
    except Exception as exc:
        raise ChemistryError(f"Meeko could not prepare the receptor: {exc}") from exc


def prepare_ligand_pdbqt(mol: Chem.Mol, seed: int) -> tuple[str, Chem.Mol]:
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        ligand = Chem.AddHs(Chem.Mol(mol))
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.useSmallRingTorsions = True
        status = AllChem.EmbedMolecule(ligand, params)
        if status != 0:
            raise ChemistryError("RDKit could not generate a ligand conformer")
        if AllChem.MMFFHasAllMoleculeParams(ligand):
            AllChem.MMFFOptimizeMolecule(ligand, maxIters=500)
        else:
            AllChem.UFFOptimizeMolecule(ligand, maxIters=500)
        setups = MoleculePreparation()(ligand)
        if not setups:
            raise ChemistryError("Meeko could not parameterize the ligand")
        pdbqt, success, message = PDBQTWriterLegacy.write_string(setups[0])
        if not success:
            raise ChemistryError(f"Meeko could not write the ligand: {message}")
        return pdbqt, ligand
    except ChemistryError:
        raise
    except Exception as exc:
        raise ChemistryError(f"ligand preparation failed: {exc}") from exc


def split_pdbqt_models(pdbqt: str, limit: int) -> list[str]:
    models = []
    current = []
    for line in pdbqt.splitlines():
        if line.startswith("MODEL"):
            current = [line]
        elif line.startswith("ENDMDL") and current:
            current.append(line)
            models.append("\n".join(current) + "\n")
            current = []
            if len(models) >= limit:
                break
        elif current:
            current.append(line)
    return models or ([pdbqt] if pdbqt.strip() else [])


def _sdf_record(mol: Chem.Mol) -> str:
    lines = [Chem.MolToMolBlock(mol).rstrip()]
    for name in mol.GetPropNames(includePrivate=False, includeComputed=False):
        lines.extend([f">  <{name}>", mol.GetProp(name), ""])
    lines.append("$$$$")
    return "\n".join(lines) + "\n"


def _float_property(mol: Chem.Mol, name: str) -> float | None:
    if not mol.HasProp(name):
        return None
    try:
        return float(mol.GetProp(name))
    except (TypeError, ValueError):
        return None


@functools.lru_cache(maxsize=1)
def gnina_version() -> str:
    binary = shutil.which(GNINA_BIN) or (GNINA_BIN if os.path.isfile(GNINA_BIN) else None)
    if not binary:
        return "unavailable"
    try:
        completed = subprocess.run(
            [binary, "--version"], capture_output=True, text=True, timeout=30, check=False,
        )
        output = (completed.stdout or completed.stderr).strip().splitlines()
        return output[0].strip() if output else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def run_gnina_docking(prepared: PreparedReceptor, prepared_ligand: Chem.Mol, request, heavy_atoms: int) -> dict:
    binary = shutil.which(GNINA_BIN) or (GNINA_BIN if os.path.isfile(GNINA_BIN) else None)
    if not binary:
        raise ChemistryError("GNINA is not installed in this chemistry-service image")

    center = prepared.box["center"]
    size = prepared.box["size"]
    with tempfile.TemporaryDirectory(prefix="terpedia-gnina-") as temporary:
        receptor_path = os.path.join(temporary, "receptor.pdbqt")
        ligand_path = os.path.join(temporary, "ligand.sdf")
        output_path = os.path.join(temporary, "poses.sdf")
        with open(receptor_path, "w", encoding="utf-8") as receptor_file:
            receptor_file.write(prepared.pdbqt)
        writer = Chem.SDWriter(ligand_path)
        writer.write(prepared_ligand)
        writer.close()

        command = [
            binary,
            "--receptor", receptor_path,
            "--ligand", ligand_path,
            "--out", output_path,
            "--center_x", str(center["x"]),
            "--center_y", str(center["y"]),
            "--center_z", str(center["z"]),
            "--size_x", str(size["x"]),
            "--size_y", str(size["y"]),
            "--size_z", str(size["z"]),
            "--scoring", request.scoring,
            "--cnn_scoring", request.cnn_scoring,
            "--pose_sort_order", "CNNscore",
            "--exhaustiveness", str(request.exhaustiveness),
            "--num_modes", str(request.poses),
            "--seed", str(request.seed),
            "--cpu", str(DOCKING_CPUS),
        ]
        if request.cnn_model == "fast":
            command.extend(["--cnn", "fast"])
        if not GNINA_USE_GPU:
            command.append("--no_gpu")
        try:
            completed = subprocess.run(
                command, capture_output=True, text=True, timeout=GNINA_TIMEOUT_SECONDS, check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise ChemistryError(f"GNINA docking exceeded {GNINA_TIMEOUT_SECONDS} seconds") from exc
        except OSError as exc:
            raise ChemistryError(f"GNINA docking could not start: {exc}") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "unknown error").strip()[-8000:]
            raise ChemistryError(f"GNINA docking failed with exit code {completed.returncode}: {detail}")
        if not os.path.exists(output_path):
            raise ChemistryError("GNINA docking completed without an output SDF")

        supplier = Chem.SDMolSupplier(output_path, removeHs=False, sanitize=False)
        output_molecules = [mol for mol in supplier if mol is not None][:request.poses]
        if not output_molecules:
            raise ChemistryError("GNINA docking produced no readable poses")

    poses = []
    for index, pose in enumerate(output_molecules):
        empirical = _float_property(pose, "minimizedAffinity")
        cnn_score = _float_property(pose, "CNNscore")
        cnn_affinity = _float_property(pose, "CNNaffinity")
        cnn_vs = _float_property(pose, "CNN_VS")
        if empirical is None or cnn_score is None:
            raise ChemistryError("GNINA output omitted required minimizedAffinity or CNNscore properties")
        poses.append({
            "rank": index + 1,
            "score_kcal_mol": round(empirical, 3),
            "cnn_score": round(cnn_score, 6),
            "cnn_affinity_pK": round(cnn_affinity, 5) if cnn_affinity is not None else None,
            "cnn_vs": round(cnn_vs, 6) if cnn_vs is not None else None,
            "ligand_efficiency_kcal_mol_per_heavy_atom": round(abs(empirical) / heavy_atoms, 4),
            "pose_sdf": _sdf_record(pose),
        })
    return {
        "poses": poses,
        "best_score_kcal_mol": poses[0]["score_kcal_mol"],
        "best_cnn_score": poses[0]["cnn_score"],
        "best_cnn_affinity_pK": poses[0]["cnn_affinity_pK"],
        "best_cnn_vs": poses[0]["cnn_vs"],
    }


def run_vina_docking(prepared: PreparedReceptor, ligand_pdbqt: str, request, heavy_atoms: int) -> dict:
    try:
        from vina import Vina

        center = list(prepared.box["center"].values())
        size = list(prepared.box["size"].values())
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdbqt") as receptor_file:
            receptor_file.write(prepared.pdbqt)
            receptor_file.flush()
            vina = Vina(sf_name=request.scoring, cpu=DOCKING_CPUS, seed=request.seed, verbosity=0)
            vina.set_receptor(receptor_file.name)
            vina.set_ligand_from_string(ligand_pdbqt)
            vina.compute_vina_maps(center=center, box_size=size)
            vina.dock(exhaustiveness=request.exhaustiveness, n_poses=max(request.poses, 10))
            energies = vina.energies(n_poses=request.poses, energy_range=request.energy_range).tolist()
            pose_models = split_pdbqt_models(
                vina.poses(n_poses=request.poses, energy_range=request.energy_range), request.poses
            )
    except Exception as exc:
        raise ChemistryError(f"AutoDock Vina docking failed: {exc}") from exc

    poses = []
    for index, energy in enumerate(energies):
        score = float(energy[0])
        poses.append({
            "rank": index + 1,
            "score_kcal_mol": round(score, 3),
            "intermolecular_kcal_mol": round(float(energy[1]), 3),
            "intramolecular_kcal_mol": round(float(energy[2]), 3),
            "torsional_kcal_mol": round(float(energy[3]), 3),
            "ligand_efficiency_kcal_mol_per_heavy_atom": round(abs(score) / heavy_atoms, 4),
            "pose_pdbqt": pose_models[index] if index < len(pose_models) else None,
        })
    return {"poses": poses, "best_score_kcal_mol": poses[0]["score_kcal_mol"] if poses else None}


async def prepare_receptor(request) -> PreparedReceptor:
    source = {"database": "user supplied"}
    preset = TARGET_PRESETS.get(request.receptor.pdb_id or "")
    reference = request.reference_ligand
    chains = list(request.receptor.chains)
    keep = list(request.receptor.keep_heterogens)

    if preset:
        if reference is None:
            from .models import ReferenceLigand
            reference = ReferenceLigand(
                residue_name=preset["reference"][0], chain=preset["reference"][1],
                residue_number=preset["reference"][2]
            )
        if not chains:
            chains = preset["chains"]
        keep = list(dict.fromkeys([*preset["keep"], *keep]))

    if request.receptor.pdb_id:
        pdb_text, source = await fetch_rcsb_receptor(request.receptor.pdb_id)
    elif request.receptor.pdb_text:
        pdb_text = request.receptor.pdb_text
    else:
        box = request.box.model_dump() if request.box else None
        return PreparedReceptor(
            request.receptor.pdbqt_text, source,
            {"engine": "user prepared PDBQT", "chains": chains, "kept_heterogens": keep}, box
        )

    if request.box:
        box = request.box.model_dump()
        box["derived_from"] = {"method": "explicit_coordinates"}
    elif reference:
        box = box_from_reference_ligand(pdb_text, reference)
    else:
        raise ChemistryError("provide a docking box or reference ligand; this PDB has no curated site preset")

    reference_mol = extract_reference_ligand(pdb_text, reference) if reference else None
    filtered, atom_count = filter_receptor_pdb(pdb_text, chains, keep, reference)
    if (preset or {}).get("pdbfixer", True):
        try:
            repaired, repair_audit = await __import__("asyncio").to_thread(
                repair_receptor_pdb, filtered, RECEPTOR_REPAIR_PH,
                int((preset or {}).get("relaxation_iterations", 20)),
            )
        except Exception as exc:
            raise ChemistryError(f"PDBFixer could not repair the receptor: {exc}") from exc
    else:
        repaired = filtered
        repair_audit = {
            "engine": "not applied",
            "reason": (preset or {}).get("repair_note"),
            "policy": "Observed protein coordinates were passed directly to Meeko; missing atoms and unresolved segments were not modeled.",
        }
    set_templates = dict((preset or {}).get("set_templates", {}))
    allow_bad_res = bool((preset or {}).get("allow_bad_res", False))
    pdbqt = await __import__("asyncio").to_thread(
        prepare_receptor_pdbqt, repaired, set_templates, allow_bad_res
    )
    return PreparedReceptor(pdbqt, source, {
        "engine": "Meeko", "engine_version": _package_version("meeko"), "rigid_receptor": True,
        "chains": chains, "kept_heterogens": keep, "prepared_atoms": atom_count,
        "explicit_residue_templates": set_templates,
        "allow_bad_res": allow_bad_res,
        "bad_residue_policy": (
            "Meeko excluded observed residues that could not match a complete chemical template."
            if allow_bad_res else "Incomplete residues are fatal."
        ),
        "removed_reference_ligand": reference.model_dump() if reference else None,
        "repair": repair_audit,
        "hydrogen_and_charge_method": f"PDBFixer hydrogens at pH {RECEPTOR_REPAIR_PH}; Meeko atom typing and charges",
        "relaxation": "Restrained local OpenMM minimization; parameters and energies are recorded in the repair audit.",
    }, box, reference_mol, reference.residue_name if reference else None)


def _package_version(name: str) -> str:
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:
        return "unknown"


async def dock_ligand(mol: Chem.Mol, label: str, request, prepared_receptor: PreparedReceptor | None = None) -> dict:
    started = time.perf_counter()
    receptor = prepared_receptor or await prepare_receptor(request)
    ligand_pdbqt, prepared_ligand = await __import__("asyncio").to_thread(prepare_ligand_pdbqt, mol, request.seed)
    heavy_atoms = prepared_ligand.GetNumHeavyAtoms()
    if request.engine == "gnina":
        docking = await __import__("asyncio").to_thread(
            run_gnina_docking, receptor, prepared_ligand, request, heavy_atoms
        )
        engine_name = "GNINA"
        engine_version = gnina_version()
        best_pose = docking["poses"][0]["pose_sdf"] if docking["poses"] else None
        pose_mime_type = "chemical/x-mdl-sdfile"
    else:
        docking = await __import__("asyncio").to_thread(
            run_vina_docking, receptor, ligand_pdbqt, request, heavy_atoms
        )
        engine_name = "AutoDock Vina"
        engine_version = _package_version("vina")
        best_pose = docking["poses"][0]["pose_pdbqt"] if docking["poses"] else None
        pose_mime_type = "chemical/x-pdbqt"
    result = {
        "ligand": {
            "label": label,
            "identifiers": identifiers(mol),
            "heavy_atoms": heavy_atoms,
            "molecular_weight": round(Descriptors.MolWt(mol), 5),
            "preparation": {"engines": ["RDKit ETKDGv3", "Meeko"], "seed": request.seed},
        },
        "receptor": {"source": receptor.source, "preparation": receptor.preparation},
        "site": receptor.box,
        "method": {
            "engine": engine_name, "engine_version": engine_version,
            "scoring_function": request.scoring, "exhaustiveness": request.exhaustiveness,
            "requested_poses": request.poses, "energy_range_kcal_mol": request.energy_range,
            **({
                "cnn_scoring": request.cnn_scoring,
                "cnn_model": request.cnn_model,
                "pose_sort_order": "CNNscore",
                "acceleration": "gpu" if GNINA_USE_GPU else "cpu",
            } if request.engine == "gnina" else {}),
        },
        **docking,
        "interpretation": {
            "score_meaning": (
                "GNINA ranks poses by CNNscore and also reports an empirical minimized affinity score."
                if request.engine == "gnina"
                else "Vina scores rank predicted poses under this exact preparation protocol."
            ),
            "not_an_affinity": "Scores are not experimental binding affinity and must not be converted to Ki, Kd, or IC50.",
            "comparison_rule": "Compare ligands only when receptor, box, protonation, scoring, and search settings are identical.",
            "validation": "Use known active/inactive controls, redocking RMSD, orthogonal scoring, and experimental binding assays.",
        },
        "sources": [
            ({
                "id": "DOCK1", "title": "GNINA molecular docking with deep learning",
                "url": "https://github.com/gnina/gnina",
            } if request.engine == "gnina" else {
                "id": "DOCK1", "title": "AutoDock Vina", "url": "https://autodock-vina.readthedocs.io/",
            }),
            {"id": "DOCK2", "title": "Meeko receptor and ligand preparation", "url": "https://meeko.readthedocs.io/"},
            {"id": "DOCK3", "title": "RCSB Protein Data Bank", "url": receptor.source.get("structure_url", "https://www.rcsb.org/")},
        ],
        "artifacts": ([{
            "type": "molecular_docking_pose", "title": f"{label} predicted pose 1",
            "mime_type": pose_mime_type, "content": best_pose,
            "receptor": receptor.source, "site": receptor.box,
        }] if best_pose else []),
        "provenance": {
            "operation": "dock", "service_version": SERVICE_VERSION,
            "rdkit_version": rdBase.rdkitVersion,
            "docking_engine": request.engine,
            "docking_engine_version": engine_version,
            "computed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "computation_ms": round((time.perf_counter() - started) * 1000, 2),
        },
    }
    return result


async def dock_ligand_batch(items: list[tuple[Chem.Mol, str, dict | None]], request) -> dict:
    started = time.perf_counter()
    receptor = await prepare_receptor(request)
    results = []
    for mol, label, resolved in items:
        result = await dock_ligand(mol, label, request, prepared_receptor=receptor)
        if resolved:
            result["ligand"]["resolution"] = resolved
        if not request.include_pose_artifacts:
            result["artifacts"] = []
            for pose in result.get("poses", []):
                pose.pop("pose_pdbqt", None)
                pose.pop("pose_sdf", None)
        results.append(result)
    return {
        "results": results,
        "count": len(results),
        "receptor": {"source": receptor.source, "preparation": receptor.preparation},
        "site": receptor.box,
        "method": {
            "engine": "GNINA" if request.engine == "gnina" else "AutoDock Vina",
            "engine_version": gnina_version() if request.engine == "gnina" else _package_version("vina"),
            "scoring_function": request.scoring, "exhaustiveness": request.exhaustiveness,
            "requested_poses": request.poses, "seed": request.seed,
            "execution": "serial ligand batch with one receptor preparation",
            **({
                "cnn_scoring": request.cnn_scoring,
                "cnn_model": request.cnn_model,
                "pose_sort_order": "CNNscore",
                "acceleration": "gpu" if GNINA_USE_GPU else "cpu",
            } if request.engine == "gnina" else {}),
        },
        "provenance": {
            "operation": "dock_batch", "service_version": SERVICE_VERSION,
            "rdkit_version": rdBase.rdkitVersion,
            "docking_engine": request.engine,
            "computed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "computation_ms": round((time.perf_counter() - started) * 1000, 2),
        },
    }


def _pose_to_rdkit(pose_text: str, pose_format: str = "pdbqt") -> Chem.Mol:
    try:
        if pose_format == "sdf":
            mol = Chem.MolFromMolBlock(pose_text.split("$$$$", 1)[0], sanitize=False, removeHs=True)
            if mol is None:
                raise ValueError("RDKit did not reconstruct the SDF pose")
            return mol
        from meeko import PDBQTMolecule, RDKitMolCreate

        pdbqt_mol = PDBQTMolecule(pose_text, is_dlg=False, skip_typing=True)
        molecules = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        mol = next((item for item in molecules if item is not None), None)
        if mol is None:
            raise ValueError("Meeko did not reconstruct a molecule")
        return Chem.RemoveHs(mol)
    except Exception as exc:
        raise ChemistryError(f"could not reconstruct docked pose for RMSD: {exc}") from exc


async def validate_cocrystal_redocking(prepared: PreparedReceptor, request) -> dict:
    if prepared.reference_ligand is None:
        return {"passed": False, "reason": "No co-crystal ligand was available for redocking."}
    result = await dock_ligand(
        prepared.reference_ligand, prepared.reference_label or "co-crystal ligand", request,
        prepared_receptor=prepared
    )
    pose_record = result.get("poses", [{}])[0]
    pose = pose_record.get("pose_sdf") or pose_record.get("pose_pdbqt")
    if not pose:
        return {"passed": False, "reason": "Redocking produced no pose."}
    docked = _pose_to_rdkit(pose, "sdf" if pose_record.get("pose_sdf") else "pdbqt")
    reference = Chem.RemoveHs(Chem.Mol(prepared.reference_ligand))
    try:
        rmsd = float(rdMolAlign.CalcRMS(reference, docked, maxMatches=100_000))
    except Exception as exc:
        raise ChemistryError(f"co-crystal RMSD calculation failed: {exc}") from exc
    threshold = float(request.redocking_rmsd_threshold)
    return {
        "passed": rmsd <= threshold, "rmsd_angstrom": round(rmsd, 4),
        "threshold_angstrom": threshold, "ligand": prepared.reference_label,
        "score_kcal_mol": result["best_score_kcal_mol"],
        **({
            "cnn_score": result.get("best_cnn_score"),
            "cnn_affinity_pK": result.get("best_cnn_affinity_pK"),
        } if request.engine == "gnina" else {}),
        "method": "symmetry-aware heavy-atom RMSD in the crystallographic coordinate frame",
    }


async def dock_hq_batch(items: list[tuple[Chem.Mol, str, dict | None]], request) -> dict:
    started = time.perf_counter()
    receptor = await prepare_receptor(request)
    validation_runs = []
    for seed in request.seeds:
        validation_request = request.model_copy(update={
            "scoring": "vina", "seed": seed,
            "exhaustiveness": request.validation_exhaustiveness,
            "poses": request.poses,
        })
        run = await validate_cocrystal_redocking(receptor, validation_request)
        validation_runs.append({"seed": seed, **run})
    passing_runs = sum(1 for run in validation_runs if run["passed"])
    required_runs = (len(validation_runs) // 2) + 1
    rmsds = [run["rmsd_angstrom"] for run in validation_runs if "rmsd_angstrom" in run]
    validation = {
        "passed": passing_runs >= required_runs,
        "passing_runs": passing_runs, "required_runs": required_runs,
        "runs": validation_runs,
        "criterion": (
            f"Top-ranked {'GNINA CNN' if request.engine == 'gnina' else 'Vina'} pose RMSD <= "
            f"{request.redocking_rmsd_threshold} A in at least {required_runs} of {len(validation_runs)} seeds."
        ),
        **({
            "median_rmsd_angstrom": round(median(rmsds), 4),
            "minimum_rmsd_angstrom": min(rmsds), "maximum_rmsd_angstrom": max(rmsds),
        } if rmsds else {}),
    }
    if request.require_redocking_pass and not validation["passed"]:
        return {
            "status": "validation_failed", "validation": validation,
            "receptor": {"source": receptor.source, "preparation": receptor.preparation},
            "site": receptor.box, "results": [],
        }

    if request.engine == "gnina":
        runs = []
        for seed in request.seeds:
            run_request = request.model_copy(update={"seed": seed})
            for mol, label, resolved in items:
                result = await dock_ligand(mol, label, run_request, prepared_receptor=receptor)
                runs.append({
                    "ligand": label,
                    "seed": seed,
                    "score_kcal_mol": result["best_score_kcal_mol"],
                    "cnn_score": result["best_cnn_score"],
                    "cnn_affinity_pK": result["best_cnn_affinity_pK"],
                    "cnn_vs": result["best_cnn_vs"],
                    "ligand_efficiency": result["poses"][0]["ligand_efficiency_kcal_mol_per_heavy_atom"],
                    "heavy_atoms": result["ligand"]["heavy_atoms"],
                    **({"resolution": resolved} if resolved else {}),
                })

        aggregates = []
        for label in sorted({run["ligand"] for run in runs}):
            ligand_runs = [run for run in runs if run["ligand"] == label]
            summary = {
                "empirical_score_kcal_mol": round(median(run["score_kcal_mol"] for run in ligand_runs), 4),
                "cnn_score": round(median(run["cnn_score"] for run in ligand_runs), 6),
                "cnn_affinity_pK": round(median(run["cnn_affinity_pK"] for run in ligand_runs), 5),
                "cnn_vs": round(median(run["cnn_vs"] for run in ligand_runs), 6),
                "ligand_efficiency_kcal_mol_per_heavy_atom": round(
                    median(run["ligand_efficiency"] for run in ligand_runs), 4
                ),
                "runs": len(ligand_runs),
            }
            aggregates.append({
                "ligand": label,
                "heavy_atoms": ligand_runs[0]["heavy_atoms"],
                "scores": summary,
            })

        ranking_fields = {
            "empirical_score_kcal_mol": False,
            "cnn_score": True,
            "cnn_affinity_pK": True,
            "cnn_vs": True,
        }
        for field, descending in ranking_fields.items():
            ranked = sorted(aggregates, key=lambda row: row["scores"][field], reverse=descending)
            for rank, row in enumerate(ranked, 1):
                row.setdefault("ranks", {})[field] = rank
        for row in aggregates:
            row["consensus_mean_rank"] = round(mean(row["ranks"].values()), 3)
        aggregates.sort(key=lambda row: (row["consensus_mean_rank"], row["ligand"]))

        return {
            "status": "completed", "validation": validation,
            "receptor": {"source": receptor.source, "preparation": receptor.preparation},
            "site": receptor.box, "results": aggregates, "runs": runs,
            "method": {
                "engine": "GNINA", "engine_version": gnina_version(),
                "cnn_scoring": request.cnn_scoring, "cnn_model": request.cnn_model,
                "seeds": list(request.seeds), "exhaustiveness": request.exhaustiveness,
                "poses": request.poses,
                "ranking": (
                    "Mean rank of median empirical energy, CNNscore, CNNaffinity, and CNN_VS. "
                    "Signal magnitudes are kept separate."
                ),
                "acceleration": "gpu" if GNINA_USE_GPU else "cpu",
            },
            "provenance": {
                "operation": "dock_hq_batch", "service_version": SERVICE_VERSION,
                "rdkit_version": rdBase.rdkitVersion, "docking_engine": "gnina",
                "docking_engine_version": gnina_version(),
                "computed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "computation_ms": round((time.perf_counter() - started) * 1000, 2),
            },
        }

    runs = []
    for scoring in request.scoring_functions:
        for seed in request.seeds:
            run_request = request.model_copy(update={"scoring": scoring, "seed": seed})
            for mol, label, resolved in items:
                result = await dock_ligand(mol, label, run_request, prepared_receptor=receptor)
                runs.append({
                    "ligand": label, "scoring": scoring, "seed": seed,
                    "score_kcal_mol": result["best_score_kcal_mol"],
                    "ligand_efficiency": result["poses"][0]["ligand_efficiency_kcal_mol_per_heavy_atom"],
                    "heavy_atoms": result["ligand"]["heavy_atoms"],
                    **({"resolution": resolved} if resolved else {}),
                })

    aggregates = []
    for label in sorted({run["ligand"] for run in runs}):
        scoring_summary = {}
        for scoring in request.scoring_functions:
            values = [run["score_kcal_mol"] for run in runs if run["ligand"] == label and run["scoring"] == scoring]
            efficiencies = [run["ligand_efficiency"] for run in runs if run["ligand"] == label and run["scoring"] == scoring]
            scoring_summary[scoring] = {
                "median_score_kcal_mol": round(median(values), 4),
                "mean_score_kcal_mol": round(mean(values), 4),
                "stdev_score_kcal_mol": round(stdev(values), 4) if len(values) > 1 else 0.0,
                "median_ligand_efficiency_kcal_mol_per_heavy_atom": round(median(efficiencies), 4),
                "best_score_kcal_mol": min(values), "runs": len(values),
            }
        ligand_runs = [run for run in runs if run["ligand"] == label]
        aggregates.append({
            "ligand": label,
            "heavy_atoms": ligand_runs[0]["heavy_atoms"],
            "scores": scoring_summary,
        })

    for scoring in request.scoring_functions:
        ranked = sorted(aggregates, key=lambda row: row["scores"][scoring]["median_score_kcal_mol"])
        for rank, row in enumerate(ranked, 1):
            row.setdefault("ranks", {})[scoring] = rank
    for row in aggregates:
        row["consensus_mean_rank"] = round(mean(row["ranks"].values()), 3)
    aggregates.sort(key=lambda row: (row["consensus_mean_rank"], row["ligand"]))

    return {
        "status": "completed", "validation": validation,
        "receptor": {"source": receptor.source, "preparation": receptor.preparation},
        "site": receptor.box, "results": aggregates, "runs": runs,
        "method": {
            "engines": list(request.scoring_functions), "seeds": list(request.seeds),
            "exhaustiveness": request.exhaustiveness, "poses": request.poses,
            "ranking": "Mean rank of per-scoring-function median scores; score magnitudes are not mixed across functions.",
        },
        "provenance": {
            "operation": "dock_hq_batch", "service_version": SERVICE_VERSION,
            "rdkit_version": rdBase.rdkitVersion,
            "computed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "computation_ms": round((time.perf_counter() - started) * 1000, 2),
        },
    }
