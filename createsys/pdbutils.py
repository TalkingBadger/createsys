from openmm.app import PDBFile
from adaptive_sampling import units


from rdkit import Chem
from rdkit.Geometry import Point3D
from ash import Fragment


from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmm.app import ForceField, PDBFile

def xyz2pdb(xyzfile: str, molname="MOL"):
    prefix = xyzfile.split('.')[0]
    frag = Fragment(xyzfile=xyzfile)
    frag.write_pdbfile_openmm(f'{prefix}.pdb', calc_connectivity=True, resname=molname)

def massFromPdb(pdbfile: str):
    return Fragment(pdbfile=pdbfile).mass

def createTopologyAndForcefieldFromPDB(pdbfile):
    #Read in coordinates of the full system
    fragment = Fragment(pdbfile=pdbfile)
    # create OpenMM system
    def add_conformer(mol, coords: list[list]):
        """Converts Bohr to Angstrom for RDKit compatibility!"""
        if len(coords) != mol.GetNumAtoms():
            raise ValueError("Number of coordinates does not match number of atoms in mol. Cannot add conformer")
        conf = Chem.Conformer(len(coords))
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(x, y, z) * units.BOHR_to_ANGSTROM)
        return mol.AddConformer(conf)

    def create_MM_openff_mol(elements: list[str], coords: list[list[float]], charge: int, bonds=None, **kwargs):
        """
        Create an RDKit Mol from a list of atomic symbols and 3D coordinates.
        
        Parameters:
            elements (List[str]): Atomic symbols, e.g., ["C", "H", "H", "H", "H]
            coords (List[List[float]]): Coordinates as [[x, y, z], ...]
            charge (int): Total charge of the molecule
            bonds (List[List[int]]): Bonds as [[0, 1], [0, 2], ...] are added as single bonds. Then bond orders are determined by RDKit.
            
        Returns:
            mol (openff.toolikt.topology.Molecule): openff molecule with embedded conformer.
        """
        from rdkit import Chem
        from rdkit.Geometry import Point3D
        from rdkit.Chem import rdDetermineBonds
        from openff.toolkit.topology import Molecule
        mol = Chem.RWMol()

        for elem in elements:
            mol.AddAtom(Chem.Atom(elem))

        add_conformer(mol, coords)

        # Add single bonds (based on topology, so bond order is lost)
        if bonds != None:
            for bond in bonds:
                mol.AddBond(bond[0], bond[1], Chem.BondType.SINGLE)
            mol = mol.GetMol()

            rdDetermineBonds.DetermineBondOrders(mol, charge, **kwargs)
        else:
            mol = mol.GetMol()
            rdDetermineBonds.DetermineBonds(mol, charge, **kwargs)
            
        return Molecule.from_rdkit(mol, hydrogens_are_explicit=True, allow_undefined_stereo=True)   

    frag = PDBFile(pdbfile)
    positions = frag.getPositions()
    topology = frag.topology
    mm_mols = []
    for chain in topology.chains():
        for res in chain.residues():
            frag = []
            elems = []
            coords = []
            for atom in res.atoms():
                frag.append(atom.index)
                elems.append(atom.element.symbol)
                pos = positions[atom.index] * 10 # nm to Angstrom
                coords.append([pos.x, pos.y, pos.z])
            if frag: # means there is an MM fragment -> Then create one
                # bonds are 0 indexing within fragment for mol generation
                bonds = [(frag.index(atom1.index), frag.index(atom2.index)) for atom1, atom2 in res.bonds()]
                charge = 0
                try:
                    mol = create_MM_openff_mol(elems, coords, charge, bonds)
                except Exception as e:
                    print("Could not create molecule for residue:", res, " with elements:", elems)
                    raise e
                mm_mols.append(mol)





    molecules = mm_mols
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules)

    # maybe more standard FFs can be added to match all possible molecules
    forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
    forcefield.registerTemplateGenerator(smirnoff.generator)


    system = forcefield.createSystem(topology)

    return topology, forcefield