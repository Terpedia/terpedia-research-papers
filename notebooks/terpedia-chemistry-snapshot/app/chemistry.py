import base64
import asyncio
import datetime as dt
import html
import json
import math
import os
from pathlib import Path
from urllib.parse import quote

import httpx
from rdkit import Chem, DataStructs, RDConfig, rdBase
from rdkit.Chem import (
    AllChem,
    ChemicalFeatures,
    Crippen,
    Descriptors,
    Draw,
    Lipinski,
    MACCSkeys,
    QED,
    rdChemReactions,
    rdDepictor,
    rdFMCS,
    rdMolDescriptors,
)
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdFingerprintGenerator


SERVICE_VERSION = "1.2.0"
MAX_ATOMS = int(os.getenv("CHEMISTRY_MAX_ATOMS", "500"))
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class ChemistryError(ValueError):
    pass


class ChemistryUpstreamError(RuntimeError):
    pass


PUBCHEM_SNAPSHOT_DATE = "2026-07-18"
CURATED_PUBCHEM_STRUCTURES = {
    "limonene": (22311, "CC1=CCC(CC1)C(=C)C"),
    "myrcene": (31253, "CC(=CCCC(=C)C=C)C"),
    "alpha-pinene": (6654, "CC1=CCC2CC1C2(C)C"),
    "beta-pinene": (14896, "CC1(C2CCC(=C)C1C2)C"),
    "linalool": (6549, "CC(=CCCC(C)(C=C)O)C"),
    "beta-caryophyllene": (5281515, "C/C/1=C\\CCC(=C)[C@H]2CC([C@@H]2CC1)(C)C"),
    "humulene": (5281520, "C/C/1=C\\CC(/C=C/C/C(=C/CC1)/C)(C)C"),
    "terpinolene": (11463, "CC1=CCC(=C(C)C)CC1"),
    "ocimene": (6434062, "CC(C)/C=C/C=C(\\C)/C=C"),
    "bisabolol": (1549992, "CC1=CC[C@@H](CC1)[C@@](C)(CCC=C(C)C)O"),
    "geraniol": (637566, "CC(=CCC/C(=C/CO)/C)C"),
    "cannabidiol": (644019, "CCCCCC1=CC(=C(C(=C1)O)[C@@H]2C=C(CC[C@H]2C(=C)C)C)O"),
    "tetrahydrocannabinol": (16078, "CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O"),
    "cannabigerol": (5315659, "CCCCCC1=CC(=C(C(=C1)O)C/C=C(\\C)/CCC=C(C)C)O"),
    "cannabinol": (2543, "CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)O"),
}
CURATED_NAME_ALIASES = {
    "pinene": "alpha-pinene",
    "caryophyllene": "beta-caryophyllene",
    "cbd": "cannabidiol",
    "thc": "tetrahydrocannabinol",
    "delta-9-thc": "tetrahydrocannabinol",
    "cbg": "cannabigerol",
    "cbn": "cannabinol",
}


def provenance(operation: str) -> dict:
    return {
        "engine": "RDKit",
        "engine_version": rdBase.rdkitVersion,
        "service_version": SERVICE_VERSION,
        "operation": operation,
        "computed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "method_note": "Computed values are model-derived molecular descriptors, not experimental measurements or clinical evidence.",
    }


def infer_format(value: str) -> str:
    stripped = value.strip()
    if stripped.startswith("InChI="):
        return "inchi"
    if "M  END" in stripped or "V2000" in stripped or "V3000" in stripped:
        return "molblock"
    return "smiles"


def parse_molecule(structure: str, structure_format: str = "auto") -> Chem.Mol:
    value = structure.strip()
    fmt = infer_format(value) if structure_format == "auto" else structure_format
    if fmt == "name":
        raise ChemistryError("name inputs must be resolved before molecule parsing")
    try:
        if fmt == "smiles":
            mol = Chem.MolFromSmiles(value)
        elif fmt == "inchi":
            mol = Chem.MolFromInchi(value)
        elif fmt == "molblock":
            mol = Chem.MolFromMolBlock(value, sanitize=True, removeHs=False)
        else:
            raise ChemistryError(f"unsupported structure format: {fmt}")
    except Exception as exc:
        raise ChemistryError(f"RDKit could not parse the {fmt} structure: {exc}") from exc
    if mol is None:
        raise ChemistryError(f"RDKit could not parse the {fmt} structure")
    if mol.GetNumAtoms() > MAX_ATOMS:
        raise ChemistryError(f"molecule has {mol.GetNumAtoms()} atoms; maximum is {MAX_ATOMS}")
    return mol


def identifiers(mol: Chem.Mol) -> dict:
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    output = {
        "canonical_smiles": canonical,
        "isomeric_smiles": canonical,
        "inchi": Chem.MolToInchi(mol),
        "inchi_key": Chem.InchiToInchiKey(Chem.MolToInchi(mol)),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
    }
    try:
        output["cxsmiles"] = Chem.MolToCXSmiles(mol)
    except Exception:
        pass
    return output


def standardization(mol: Chem.Mol) -> dict:
    cleanup = rdMolStandardize.Cleanup(mol)
    parent = rdMolStandardize.FragmentParent(cleanup)
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent)
    enumerator = rdMolStandardize.TautomerEnumerator()
    tautomer = enumerator.Canonicalize(uncharged)
    return {
        "clean_smiles": Chem.MolToSmiles(cleanup, isomericSmiles=True),
        "fragment_parent_smiles": Chem.MolToSmiles(parent, isomericSmiles=True),
        "uncharged_parent_smiles": Chem.MolToSmiles(uncharged, isomericSmiles=True),
        "canonical_tautomer_smiles": Chem.MolToSmiles(tautomer, isomericSmiles=True),
    }


def descriptors(mol: Chem.Mol) -> dict:
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotors = Lipinski.NumRotatableBonds(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    return {
        "molecular_weight": round(mw, 5),
        "exact_mass": round(Descriptors.ExactMolWt(mol), 5),
        "logp": round(logp, 5),
        "tpsa": round(tpsa, 5),
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": rotors,
        "heavy_atoms": Lipinski.HeavyAtomCount(mol),
        "ring_count": Lipinski.RingCount(mol),
        "aromatic_rings": Lipinski.NumAromaticRings(mol),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 5),
        "formal_charge": formal_charge,
        "radical_electrons": Descriptors.NumRadicalElectrons(mol),
        "stereocenters": [{"atom": atom, "assignment": assignment} for atom, assignment in chiral],
        "bertz_complexity": round(Descriptors.BertzCT(mol), 5),
        "qed": round(QED.qed(mol), 5),
        "rules": {
            "lipinski_rule_of_five": {
                "passes": sum([mw > 500, logp > 5, hbd > 5, hba > 10]) <= 1,
                "violations": sum([mw > 500, logp > 5, hbd > 5, hba > 10]),
            },
            "veber": {"passes": rotors <= 10 and tpsa <= 140, "rotatable_bonds_max": 10, "tpsa_max": 140},
            "rule_of_three": {"passes": mw <= 300 and logp <= 3 and hbd <= 3 and hba <= 3 and rotors <= 3 and tpsa <= 60},
        },
    }


def alert_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def structural_alerts(mol: Chem.Mol) -> list[dict]:
    matches = alert_catalog().GetMatches(mol)
    unique = {}
    for match in matches:
        description = match.GetDescription()
        unique[description] = {"description": description}
    return list(unique.values())


def molecule_svg(mol: Chem.Mol, legend: str = "", width: int = 600, height: int = 420,
                 highlight_atoms: list[int] | None = None) -> str:
    drawable = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(drawable)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().legendFontSize = 18
    drawer.DrawMolecule(drawable, legend=legend, highlightAtoms=highlight_atoms or [])
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace("svg:", "")


def svg_artifact(svg: str, title: str, artifact_type: str = "molecule_2d") -> dict:
    return {
        "type": artifact_type,
        "title": title,
        "mime_type": "image/svg+xml",
        "encoding": "utf-8",
        "content": svg,
        "alt": f"RDKit {artifact_type.replace('_', ' ')} illustration for {title}",
    }


def generate_conformers(mol: Chem.Mol, count: int) -> dict | None:
    if count <= 0:
        return None
    work = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF8A1
    params.pruneRmsThresh = 0.35
    params.numThreads = 0
    conformer_ids = list(AllChem.EmbedMultipleConfs(work, numConfs=count, params=params))
    if not conformer_ids:
        raise ChemistryError("RDKit could not generate a 3D conformer")
    energies = []
    method = "UFF"
    if AllChem.MMFFHasAllMoleculeParams(work):
        method = "MMFF94"
        results = AllChem.MMFFOptimizeMoleculeConfs(work, numThreads=0, maxIters=500)
    else:
        results = AllChem.UFFOptimizeMoleculeConfs(work, numThreads=0, maxIters=500)
    for conformer_id, (status, energy) in zip(conformer_ids, results):
        energies.append({"conformer_id": conformer_id, "status": status, "energy_kcal_mol": round(float(energy), 6)})
    energies.sort(key=lambda item: item["energy_kcal_mol"])
    best_id = energies[0]["conformer_id"]
    molblock = Chem.MolToMolBlock(work, confId=best_id)
    return {
        "method": f"ETKDGv3 + {method}",
        "requested": count,
        "generated": len(conformer_ids),
        "lowest_energy_conformer_id": best_id,
        "relative_energies_kcal_mol": [
            {**entry, "relative_energy_kcal_mol": round(entry["energy_kcal_mol"] - energies[0]["energy_kcal_mol"], 6)}
            for entry in energies
        ],
        "sdf_molblock": molblock,
        "sdf_base64": base64.b64encode(molblock.encode()).decode(),
        "caveat": "Force-field energies rank generated conformers only; they are not binding free energies or solution-state populations.",
    }


def analyze_molecule(mol: Chem.Mol, label: str = "Molecule", depict: bool = True,
                     conformers: int = 0, include_alerts: bool = True) -> dict:
    ids = identifiers(mol)
    result = {
        "label": label,
        "identifiers": ids,
        "standardization": standardization(mol),
        "descriptors": descriptors(mol),
        "structural_alerts": structural_alerts(mol) if include_alerts else [],
        "provenance": provenance("analyze"),
        "artifacts": [],
    }
    if depict:
        result["artifacts"].append(svg_artifact(molecule_svg(mol, label), label))
    conformer_result = generate_conformers(mol, conformers)
    if conformer_result:
        result["conformers"] = conformer_result
    return result


def fingerprint(mol: Chem.Mol, kind: str, radius: int, n_bits: int):
    if kind == "morgan":
        return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits).GetFingerprint(mol)
    if kind == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    return Chem.RDKFingerprint(mol, fpSize=n_bits)


def compare_molecules(items: list[tuple[str, Chem.Mol]], kind: str, radius: int,
                      n_bits: int, include_mcs: bool, depict: bool) -> dict:
    fps = [fingerprint(mol, kind, radius, n_bits) for _, mol in items]
    matrix = []
    pairs = []
    for i, (label, _) in enumerate(items):
        row = []
        for j in range(len(items)):
            score = float(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
            row.append(round(score, 6))
            if i < j:
                pairs.append({"a": items[i][0], "b": items[j][0], "tanimoto": round(score, 6)})
        matrix.append(row)
    result = {
        "labels": [label for label, _ in items],
        "fingerprint": {"type": kind, "radius": radius if kind == "morgan" else None, "bits": n_bits},
        "similarity_matrix": matrix,
        "pairs": sorted(pairs, key=lambda pair: pair["tanimoto"], reverse=True),
        "provenance": provenance("compare"),
        "artifacts": [],
    }
    if include_mcs:
        mcs = rdFMCS.FindMCS([mol for _, mol in items], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
        result["maximum_common_substructure"] = {
            "smarts": mcs.smartsString,
            "atoms": mcs.numAtoms,
            "bonds": mcs.numBonds,
            "canceled": mcs.canceled,
        }
    if depict:
        image = Draw.MolsToGridImage(
            [mol for _, mol in items], molsPerRow=min(4, len(items)), subImgSize=(300, 240),
            legends=[label for label, _ in items], useSVG=True
        )
        result["artifacts"].append(svg_artifact(str(image), "Molecular comparison", "molecule_grid"))
    return result


def substructure_search(query: Chem.Mol, targets: list[tuple[str, Chem.Mol]], use_chirality: bool,
                        depict: bool) -> dict:
    matches = []
    artifacts = []
    for label, target in targets:
        atom_matches = target.GetSubstructMatches(query, useChirality=use_chirality, maxMatches=1000)
        matches.append({"label": label, "matched": bool(atom_matches), "match_count": len(atom_matches), "atom_matches": atom_matches})
        if depict and atom_matches:
            artifacts.append(svg_artifact(molecule_svg(target, label, highlight_atoms=list(atom_matches[0])), label, "substructure_match"))
    return {
        "query_smarts": Chem.MolToSmarts(query),
        "use_chirality": use_chirality,
        "matches": matches,
        "artifacts": artifacts,
        "provenance": provenance("substructure"),
    }


def reaction_analysis(value: str, reaction_format: str, depict: bool) -> dict:
    try:
        reaction = rdChemReactions.ReactionFromSmarts(value, useSmiles=reaction_format == "smiles")
    except Exception as exc:
        raise ChemistryError(f"RDKit could not parse the reaction: {exc}") from exc
    if reaction is None:
        raise ChemistryError("RDKit could not parse the reaction")
    result = {
        "reactants": reaction.GetNumReactantTemplates(),
        "products": reaction.GetNumProductTemplates(),
        "agents": reaction.GetNumAgentTemplates(),
        "canonical_smarts": rdChemReactions.ReactionToSmarts(reaction),
        "canonical_smiles": rdChemReactions.ReactionToSmiles(reaction),
        "artifacts": [],
        "provenance": provenance("reaction"),
        "caveat": "A parsed reaction scheme does not predict yield, selectivity, feasibility, or safety.",
    }
    if depict:
        drawer = rdMolDraw2D.MolDraw2DSVG(1000, 360)
        drawer.DrawReaction(reaction)
        drawer.FinishDrawing()
        result["artifacts"].append(svg_artifact(drawer.GetDrawingText(), "Reaction scheme", "reaction_scheme"))
    return result


def pharmacophore_features(mol: Chem.Mol, count: int) -> dict:
    conformer_data = generate_conformers(mol, count)
    work = Chem.MolFromMolBlock(conformer_data["sdf_molblock"], removeHs=False)
    feature_factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    features = []
    for feature in feature_factory.GetFeaturesForMol(work):
        position = feature.GetPos()
        features.append({
            "family": feature.GetFamily(),
            "type": feature.GetType(),
            "atom_ids": list(feature.GetAtomIds()),
            "position": {"x": round(position.x, 5), "y": round(position.y, 5), "z": round(position.z, 5)},
        })
    return {
        "features": features,
        "conformer": conformer_data,
        "provenance": provenance("pharmacophore"),
        "caveat": "Feature points are ligand-derived hypotheses; they do not establish target binding or biological activity.",
    }


async def resolve_pubchem_name(name: str) -> dict:
    normalized = " ".join(name.strip().lower().split())
    canonical_name = CURATED_NAME_ALIASES.get(normalized, normalized)
    curated = CURATED_PUBCHEM_STRUCTURES.get(canonical_name)
    if curated:
        cid, smiles = curated
        mol = parse_molecule(smiles, "smiles")
        ids = identifiers(mol)
        source_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        return {
            "query": name,
            "cid": cid,
            "title": canonical_name,
            "smiles": smiles,
            "inchi": ids["inchi"],
            "inchi_key": ids["inchi_key"],
            "formula": ids["formula"],
            "molecular_weight": str(round(Descriptors.MolWt(mol), 5)),
            "resolution_method": "curated_pubchem_snapshot",
            "snapshot_date": PUBCHEM_SNAPSHOT_DATE,
            "source": {"title": "PubChem", "url": source_url, "retrieved_at": PUBCHEM_SNAPSHOT_DATE},
        }
    properties = "Title,CanonicalSMILES,IsomericSMILES,ConnectivitySMILES,SMILES,InChI,InChIKey,MolecularFormula,MolecularWeight"
    source_url = f"{PUBCHEM_BASE}/compound/name/{quote(name, safe='')}/property/{properties}/JSON"
    response = None
    async with httpx.AsyncClient(timeout=12, follow_redirects=True) as client:
        for attempt in range(3):
            response = await client.get(source_url, headers={"User-Agent": "Terpedia-Chemistry/1.0"})
            if response.status_code not in {429, 500, 502, 503, 504}:
                break
            if attempt < 2:
                await asyncio.sleep(0.4 * (2 ** attempt))
    if response.status_code == 404:
        raise ChemistryError(f"PubChem did not resolve {name!r}")
    if response.status_code in {429, 500, 502, 503, 504}:
        raise ChemistryUpstreamError(f"PubChem is temporarily unavailable ({response.status_code})")
    response.raise_for_status()
    rows = response.json().get("PropertyTable", {}).get("Properties", [])
    if not rows:
        raise ChemistryError(f"PubChem returned no structure for {name!r}")
    row = rows[0]
    smiles = row.get("IsomericSMILES") or row.get("SMILES") or row.get("CanonicalSMILES") or row.get("ConnectivitySMILES")
    if not smiles:
        raise ChemistryError(f"PubChem returned no SMILES for {name!r}")
    return {
        "query": name,
        "cid": row.get("CID"),
        "title": row.get("Title") or name,
        "smiles": smiles,
        "inchi": row.get("InChI"),
        "inchi_key": row.get("InChIKey"),
        "formula": row.get("MolecularFormula"),
        "molecular_weight": row.get("MolecularWeight"),
        "resolution_method": "pubchem_live",
        "source": {"title": "PubChem", "url": source_url, "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat()},
    }


def compact_for_agent(result: dict) -> str:
    without_artifacts = {key: value for key, value in result.items() if key != "artifacts"}
    return json.dumps(without_artifacts, separators=(",", ":"), default=str)
