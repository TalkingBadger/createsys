import os
import numpy as np
from . import packmolBuilder, packmolStructure
from .pdbutils import createTopologyAndForcefieldFromPDB
import ash
from adaptive_sampling import units

def calculate_sphere_radius(mass, density=1):
    """Calculates the radius of a sphere based on mass and density (Angstrom).
    Args:
        mass (float): mass in amu
        density (float): density in g/cm^3
    Returns:
        radius (float): radius in Angstrom
    """
    mass_g = mass * units.atomic_to_kg * 1000  # amu -> g
    V = mass_g / density
    return (V / np.pi * 0.75) ** (1 / 3) * 1e8  # cm -> Angstrom

def calculate_n_solv(mass_solv, boxsize, radius_inner_sphere, density=1, rounded=True):
    """Calculate number of solvent molecules to match given density for the volume around the inner sphere.
    Args:
        mass_solv (float): mass of the current solvent molecules in amu
        boxsize (float): cubic box side length in Angstrom
        radius_inner_sphere (float): radius of the inner sphere in Angstrom
        density (float): target density in g/cm^3
        rounded (bool): if True, return rounded integer number of molecules
    Returns:
        n_solv (int or float): multiple of current solvent molecules needed to match density
    """
    V_sphere = 4/3 * np.pi * radius_inner_sphere ** 3
    V = boxsize ** 3 - V_sphere
    mass_g = V * density * 1e-24  # A^3 -> cm^3
    n_solv = mass_g / (mass_solv * units.atomic_to_kg * 1000)
    if rounded:
        return int(np.round(n_solv))
    return n_solv

def calculate_mass(mols: dict):
    """Calculate total mass of given molecules dictionary.

    Args:
        mols (dict): dictionary of molecule names and their numbers.
    """
    return sum(
        num * ash.Fragment(pdbfile=f"{mol}.pdb", printlevel=0).mass
        for mol, num in mols.items()
    )
    
def scale_n_solv(mols_solvent: dict, boxsize: float, radius_inner_sphere: float, density: float = 1):
    """Scale number of solvent molecules to match density in the box.
    Args:
        mols_solvent (dict): dictionary of solvent molecule names and their numbers.
        boxsize (float): cubic box side length in Angstrom
        radius_inner_sphere (float): radius of the inner sphere in Angstrom
        density (float): target density in g/cm^3
    Returns:
        mols_solvent_new (dict): new dictionary of solvent molecule names and their scaled numbers.
    """
    mass = calculate_mass(mols_solvent)
    n_solv = calculate_n_solv(mass, boxsize, radius_inner_sphere=radius_inner_sphere, density=density)
    mols_solvent_new = {} 
    for mol, num in mols_solvent.items():
        mols_solvent_new[mol] = int(num * n_solv)
    return mols_solvent_new

def rename_residues_inner_outer(outfile: str, num_inner_species: int, inner_name: str = "INN", outer_name: str = "OUT"):
    """Rename residues in PDB file for easier identification of inner and outer molecules.
    Args:
        outfile (str): output PDB filename
        num_inner_species (int): number of different inner molecule species
        inner_name (str): residue name for inner molecules
        outer_name (str): residue name for outer molecules
    Returns:
        None
    """
    import mdtraj as md
    traj = md.load_pdb(outfile)
    inner_name = "INN"
    outer_name = "OUT"
    for residue in traj.topology.residues:
        if residue.chain.index < num_inner_species:
            residue.name = inner_name
        else:
            residue.name = outer_name
    traj.save_pdb(outfile)

def generate_initial_config(inner_mols: dict, outer_mols: dict, seed: int, radius_sphere: float, boxsize: float = 30, solvent_outside: bool = True, outfile: str = "system.pdb", rename_residues: bool = True):
    """
    Generate initial configuration with Packmol for arbitrary inner and outer molecules.

    Args:
        inner_mols (dict): molecules inside the inner sphere, e.g. {"adenine": 5, "h2o": 10}
        outer_mols (dict): molecules in the outer environment, e.g. {"h2o": 50}
        seed (int): random seed for Packmol
        radius_sphere (float): inner sphere radius in Angstrom
        boxsize (float): cubic box side length in Angstrom
        solvent_outside (bool): if True, outer molecules will only be placed outside the inner sphere
        outfile (str): output PDB filename
        rename_residues (bool): if True, rename residues for easier identification into INN and OUT depending on region
    Returns:
        None
    """

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

    # Packmol builder
    pb = packmolBuilder(outfile=outfile)
    pb.set_seed(seed)
    pb.set_pbc(boxsize)

    # Add inner molecules
    for mol, num in inner_mols.items():
        m = packmolStructure(f"{mol}.pdb", num)
        m.inside_sphere(radius_sphere, boxsize / 2)
        pb.add_structure(m)

    # Add outer molecules
    for mol, num in outer_mols.items():
        m = packmolStructure(f"{mol}.pdb", num)
        m.inside_box(boxsize)
        if solvent_outside:
            m.outside_sphere(radius_sphere, boxsize / 2)
        pb.add_structure(m)

    # Create input file and run CLI packmol
    pb.run(outfile)

    if rename_residues:
        rename_residues_inner_outer(outfile, num_inner_species=len(inner_mols))
        
