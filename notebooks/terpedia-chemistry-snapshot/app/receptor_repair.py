import io
from importlib.metadata import version


def _atom_names(atoms):
    return [atom.name if hasattr(atom, "name") else str(atom) for atom in atoms]


def _observed_atom_keys(lines):
    return {
        (line[21:22].strip(), line[22:26].strip(), line[12:16].strip())
        for line in lines if len(line) >= 26
    }


def _relax_rebuilt_atoms(fixer, observed_atom_keys, max_iterations):
    from openmm import CustomExternalForce, LocalEnergyMinimizer, VerletIntegrator
    from openmm.app import CutoffNonPeriodic, ForceField, HBonds, Simulation
    from openmm.unit import kilojoule_per_mole, nanometer, picoseconds

    forcefield = ForceField("amber14-all.xml")
    system = forcefield.createSystem(
        fixer.topology, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.0 * nanometer, constraints=HBonds,
    )
    restraint = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 1000.0)
    for parameter in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(parameter)
    restrained = 0
    for atom in fixer.topology.atoms():
        key = (atom.residue.chain.id, atom.residue.id, atom.name)
        if key not in observed_atom_keys or atom.element.symbol == "H":
            continue
        position = fixer.positions[atom.index].value_in_unit(nanometer)
        restraint.addParticle(atom.index, [position.x, position.y, position.z])
        restrained += 1
    system.addForce(restraint)
    integrator = VerletIntegrator(0.001 * picoseconds)
    simulation = Simulation(fixer.topology, system, integrator)
    simulation.context.setPositions(fixer.positions)
    initial = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    LocalEnergyMinimizer.minimize(simulation.context, tolerance=50.0, maxIterations=max_iterations)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    fixer.positions = state.getPositions()
    final = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    return {
        "method": "restrained OpenMM minimization", "restrained_observed_heavy_atoms": restrained,
        "max_iterations": max_iterations, "initial_potential_kj_mol": round(initial, 3),
        "final_potential_kj_mol": round(final, 3),
    }


def repair_receptor_pdb(pdb_text: str, ph: float = 7.4, relaxation_iterations: int = 20) -> tuple[str, dict]:
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    heterogens = [line for line in pdb_text.splitlines() if line.startswith("HETATM")]
    # Filter-generated streams can contain an orphan TER after selecting a
    # chain or removing a ligand. Emit canonical protein coordinates for
    # OpenMM; chain IDs remain encoded in each ATOM record.
    protein_lines = [line for line in pdb_text.splitlines() if line.startswith("ATOM  ")]
    if not protein_lines:
        raise ValueError("receptor contains no protein ATOM records")
    observed_atom_keys = _observed_atom_keys(protein_lines)
    protein_lines.append("END")
    fixer = PDBFixer(pdbfile=io.StringIO("\n".join(protein_lines) + "\n"))
    fixer.findMissingResidues()
    missing_residues = [
        {"chain_index": chain, "residue_index": index, "residue_names": list(names)}
        for (chain, index), names in sorted(fixer.missingResidues.items())
    ]
    # Loop building can change a binding pocket substantially. Record gaps but
    # do not model unresolved residues in this rigid-receptor workflow.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstandard = [
        {"chain": residue.chain.id, "residue": residue.id, "name": residue.name, "replacement": replacement}
        for residue, replacement in fixer.nonstandardResidues
    ]
    if nonstandard:
        fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    missing_atoms = [
        {"chain": residue.chain.id, "residue": residue.id, "name": residue.name, "atoms": _atom_names(atoms)}
        for residue, atoms in fixer.missingAtoms.items()
    ]
    missing_terminals = [
        {"chain": residue.chain.id, "residue": residue.id, "name": residue.name, "atoms": _atom_names(atoms)}
        for residue, atoms in fixer.missingTerminals.items()
    ]
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    try:
        relaxation = _relax_rebuilt_atoms(fixer, observed_atom_keys, relaxation_iterations)
    except Exception as exc:
        relaxation = {"method": "not applied", "error": str(exc)}
    output = io.StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, output, keepIds=True)
    repaired_lines = [line for line in output.getvalue().splitlines() if line != "END"]
    serials = []
    for line in repaired_lines:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                serials.append(int(line[6:11]))
            except ValueError:
                pass
    next_serial = max(serials, default=0) + 1
    preserved_heterogens = []
    for line in heterogens:
        preserved_heterogens.append(f"{line[:6]}{next_serial:>5}{line[11:]}")
        next_serial += 1
    repaired_lines.extend(preserved_heterogens)
    repaired_lines.append("END")
    return "\n".join(repaired_lines) + "\n", {
        "engine": "PDBFixer", "engine_version": version("pdbfixer"), "ph": ph,
        "missing_residue_segments_not_built": missing_residues,
        "replaced_nonstandard_residues": nonstandard,
        "added_missing_heavy_atoms": missing_atoms,
        "added_terminal_atoms": missing_terminals,
        "preserved_heterogen_atoms": len(preserved_heterogens),
        "heterogen_policy": "Heterogens were excluded from PDBFixer and restored unchanged after protein repair.",
        "relaxation": relaxation,
        "policy": "Missing atoms and hydrogens were added; unresolved residue segments were recorded but not modeled.",
    }
