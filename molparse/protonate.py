import mrich as logger
from pathlib import Path

def protonate(
    sys: "System",
    minimise: bool = False,
    pH: float = 7.0,
    remove_residues: list[str] | None = ["DMS", "TRS", "LIG", "CL"],
    trim_terminal_residues: int = 0,
    out_file: str | Path = None,
    return_file: bool = False,
) -> "System | str":
    """Protonate a system using pdbfixer and openmm

    :param sys: Input :class:`.System`
    :param minimise: Perform an energy minimisation?
    :param pH: System pH
    :param remove_residues: list of residue names to remove
    :param out_file: write the protonated PDB here
    :param trim_terminal_residues: number of residues to trim from each end of all protein chains
    :returns: Output :class:`.System`
    """

    from .amber import prep4amber
    from .io import parsePDB, writePDB
    from openmm.app import PME, ForceField, Modeller, PDBFile, Simulation
    from tempfile import NamedTemporaryFile
    from pdbfixer import PDBFixer

    orig_sys = sys.copy()

    if remove_residues:
        orig_sys.remove_residues(names=remove_residues, no_summary=True)

    if trim_terminal_residues:
        for chain in orig_sys.chains:
            if chain.type != "PRO":
                continue
            residues = [r for r in chain.residues[:trim_terminal_residues]] + [
                r for r in chain.residues[-trim_terminal_residues:]
            ]
            indices = [r.index for r in residues]
            chain.remove_residues(indices=indices)

    # prepare IO files
    pdb_orig = NamedTemporaryFile(mode="w+t", suffix=".pdb")
    pdb_prot = NamedTemporaryFile(mode="w+t", suffix=".pdb")

    # write the original
    writePDB(pdb_orig.name, orig_sys, verbosity=0)
    # writePDB("orig.pdb", orig_sys, verbosity=0)

    # fix the PDB termini and hydrogens
    fixer = PDBFixer(pdb_orig.name)
    fixer.findMissingResidues()
    if fixer.missingResidues:
        logger.warning(f"{fixer.missingResidues=}")
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        logger.warning(f"{fixer.nonstandardResidues=}")
    fixer.findMissingAtoms()
    if fixer.missingAtoms:
        logger.warning(f"{fixer.missingAtoms=}")
    if fixer.missingTerminals:
        logger.print(f"{fixer.missingTerminals=}")
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(pH)
    PDBFile.writeFile(fixer.topology, fixer.positions, pdb_prot)
    # PDBFile.writeFile(fixer.topology, fixer.positions, "prot.pdb")

    if minimise:

        # openmm implementation is sensitive...

        pdb_mini = NamedTemporaryFile(mode="w+t", suffix=".pdb")

        # from rdkit import Chem
        # from rdkit.Chem import AllChem
        
        # mol = Chem.MolFromPDBFile(pdb_prot.name, removeHs=False)
        # # Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
        # # mol = Chem.AddHs(mol, addCoords=False)
        # # mol.UpdatePropertyCache(strict=False)
        
        # # use 'MMFF94s' for biopolymers; fall back to UFF if unavailable
        # try:
        #     mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        # except Exception as e:
        #     logger.warning(e)
        #     mp = None
        #     raise
            
        # if mp is not None:
        #     ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=0)
        # else:
        #     logger.warning("MMFF not available. Falling back to UFF")
        #     ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        
        # for atom in mol.GetAtoms():
        #     if atom.GetAtomicNum() != 1:   # fix all heavy atoms
        #         ff.AddFixedPoint(atom.GetIdx())
        
        # ff.Minimize(maxIts=1_000)
        # Chem.MolToPDBFile(mol, pdb_mini.name)
        # Chem.MolToPDBFile(mol, "mini.pdb")

        # openmm objects
        pdb = PDBFile(pdb_prot.name)
        
        # forcefield = ForceField("amber99sb.xml", "tip3p.xml")
        forcefield = ForceField("amber14-all.xml")

        modeller = Modeller(pdb.topology, pdb.positions)

        modeller.addHydrogens(pH=pH)

        # prepare the simulation
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        integrator = VerletIntegrator(0.001 * picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        # minimise proton positions
        simulation.minimizeEnergy(maxIterations=100)

        # update positions
        positions = simulation.context.getState(getPositions=True).getPositions()

        # write files
        PDBFile.writeFile(simulation.topology, positions, pdb_mini)

        # read in mp.System
        sys = parsePDB(pdb_mini.name, verbosity=0)
        pdb_mini.close()

    else:
        # read in mp.System
        sys = parsePDB(pdb_prot.name, verbosity=0)
        pdb_prot.close()

    # fix residue numbering and chain naming
    for orig_chain, new_chain in zip(orig_sys.chains, sys.chains):

        if orig_chain.sequence != new_chain.sequence:

            assert len(orig_chain.sequence) > len(
                new_chain.sequence
            ), f"Sequences don't match: {orig_chain.sequence} {new_chain.sequence}"

            offset = orig_chain.sequence.find(new_chain.sequence)
            assert offset >= 0

        else:
            offset = 0

        new_chain.name = orig_chain.name

        for i, new_residue in enumerate(new_chain.residues):
            orig_residue = orig_chain.residues[i + offset]
            assert orig_residue.name == new_residue.name
            new_residue.number = orig_residue.number

    # close files
    pdb_orig.close()
    
    if out_file:
        pdb_out = Path(out_file)
        assert pdb_out.name.endswith(".pdb")
        logger.writing(pdb_out)
        sys.write(pdb_out, verbosity=0)
        
    if return_file:
        if out_file:
            return sys, pdb_out
        else:
            pdb_out = NamedTemporaryFile(mode="w+t", suffix=".pdb").name
            sys.write(pdb_out.name, verbosity=0)
            return sys, pdb_out
        
    return sys
