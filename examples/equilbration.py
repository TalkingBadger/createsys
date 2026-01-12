from createsys import Equilibrator


# now run the initial_configuration example to create a system
from initial_configuration import outfile, radius_inner

eq = Equilibrator(outfile)

eq.minimize_energy(max_iterations=500)
eq.warmup_simulation()

eq.equilibrate_box(numsteps_per_NPT=1000, volume_threshold=2, density_threshold=0.1)

eq.save_positions('equilibrated_system.pdb')