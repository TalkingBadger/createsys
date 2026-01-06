from createsys.initial_config_utils import calculate_sphere_radius, calculate_mass, scale_n_solv, generate_initial_config

# This is an example script to demonstrate the usage of the initial configuration utilities.

# Often, you know the exact number of molecules you want inside the inner sphere, 
# but the box should be filled with solvent to a specific density.
# If you know the exact number of molecules and the desired radius of the inner sphere, you can directly provide them to `generate_initial_config`.

# The inner molecules. Molecule names should correspond to available PDB files (mol.pdb).
mols_inner = {
    "adenine": 5,
    "acetaldehyde": 5,
    "d_glyceraldehyde": 5,
    "h2o": 15,
}

# The solvent molecules. Usually, water is enough, but you can add solvent mixtures and scale them accordingly.
mols_solvent = {
    "h2o": 2,
    "acetaldehyde": 1,
}

# Define targeted density (g/cm^3), and cubical box size (Angstrom).
target_density = 1.0  # g/cm^3
box_size = 30.0  # Angstrom

# Now calculate the radius of the inner sphere to match density.
mass_inner = calculate_mass(mols_inner)
radius_inner = calculate_sphere_radius(mass_inner, density=target_density)

# Scale the number of solvent molecules to match the target density. The ratio will be preserved.
mols_solvent_scaled = scale_n_solv(mols_solvent, boxsize=box_size, radius_inner_sphere=radius_inner, density=target_density)

# Now generate the initial configuration with the calculated parameters. Packmol will be used under the hood.
outfile = "example_system_adenine.pdb"
generate_initial_config(
    inner_mols=mols_inner,
    outer_mols=mols_solvent_scaled,
    seed=1,
    radius_sphere=radius_inner,
    boxsize=box_size,
    solvent_outside=True,           # Solvent will only be added arround the inner sphere
    outfile="example_system_adenine.pdb",
    rename_residues=True,         # Rename residues for easier identification into INN and OUT    
)   

# Note that you can use the chain IDs to identify the different molecular species 
# and the resnames to identify the regions.
