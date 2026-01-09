import os
import numpy as np
import ash
from openmm import app, unit
from adaptive_sampling.interface import AdaptiveSamplingOpenMM
from openmm.unit import Quantity
from openmm import MonteCarloBarostat

from adaptive_sampling.interface import AdaptiveSamplingOpenMM
from adaptive_sampling.exploration_tools import QMMMHyperreactor

from .pdbutils import createTopologyAndForcefieldFromPDB, getAtomIndicesFromResidueNames


class Equilibrator:
    """
        Class to handle equilibration of a system created with createsys. Equilibration is done with `openmm` and `AdaptiveSamplingOpenMM`.
    """
    def __init__(self, pdbfile, default_temperature=300.0, default_dt = 2.0, save_traj=False):
        """
        Args:
            pdbfile (str): input PDB file with initial configuration
            default_temperature (float): default temperature for MD in K
            default_dt (float): default time step for MD in fs
            save_traj (bool): if True, save trajectory files during equilibration
        """
        self.topology, self.forcefield = createTopologyAndForcefieldFromPDB(pdbfile)
        self.indices_inner = getAtomIndicesFromResidueNames(pdbfile, ['INN'])
        self.pdb = app.PDBFile(pdbfile) 
        self.positions = self.pdb.positions
        # Handle periodic box correctly
        box_vectors = self.pdb.topology.getPeriodicBoxVectors()
        self.topology.setPeriodicBoxVectors(box_vectors)
        self.system = self.forcefield.createSystem(
            self.topology,
            constraints=app.HBonds,   # constrain all bonds to hydrogens
            nonbondedMethod=app.PME,          # Particle Mesh Ewald (periodic)
            nonbondedCutoff=1.2*unit.nanometer,
        )
        self.default_temperature = default_temperature
        self.default_dt = default_dt
        self.save_traj = save_traj
        self._inner_constraint = {
            'relative_size': 0.99,
            'k_conf': 0.0,
        }
      
    def set_inner_constraint(self, relative_size: float, k_conf: float = 0.25):
        """Add spherical constraint to inner region to keep it together during equilibration.
        
        Args:
            relative_size (float): relative size of the spherical constraint compared to half box size
            k_conf (float): force constant of the spherical constraint in kJ/mol/A^2
        """
        self._inner_constraint = {
            'relative_size': relative_size,
            'k_conf': k_conf
        } 
        
    def add_spherical_constraint_to_inner(self, md: AdaptiveSamplingOpenMM, relative_size: float, k_conf: float = 0.25):
        """Add spherical constraint to inner region to keep it together during equilibration.
        Acts directly on an instancel of AdaptiveSamplingOpenMM.
        Args:
            md (AdaptiveSamplingOpenMM): MD object to which the constraint will be added
            relative_size (float): relative size of the spherical constraint potential free radius compared to half box size
            k_conf (float): force constant of the spherical constraint in kJ/mol/A^2
        """ 
        md.qmatoms = self.indices_inner

        bias_pot = QMMMHyperreactor(
            k_conf=k_conf,
            relative_soft_box_size=relative_size,
            qm_confinement='spherical',
            override_qm_indices=self.indices_inner,
            amd_parameter=0,
            init_steps=1e10,
            equil_steps=1e10,
            md=md
        )
        md.set_sampling_algorithm(bias_pot)
    
    def create_MD(self, dt=2.0, equil_temp=300.0, langevin_damping=1.0, cv_atoms=[], **kwargs):
        """Create AdaptiveSamplingOpenMM MD object for equilibration.\
        Args:
            dt (float): MD timestep in fs
            equil_temp (float): equilibrium temperature for langevin dynamics in K
            langevin_damping (float): friction factor for langevin dynamics in 1/ps
            cv_atoms (list): indices of atoms that are involved in Collective Variable
            **kwargs: additional arguments for AdaptiveSamplingOpenMM
        Returns:
            md (AdaptiveSamplingOpenMM): MD object for equilibration
        """
        md = AdaptiveSamplingOpenMM(
            positions=self.positions,
            topology=self.topology,
            system=self.system,
            dt=dt,
            cv_atoms=cv_atoms,
            equil_temp=equil_temp,
            langevin_damping=langevin_damping,
            **kwargs,
        )
        self.add_spherical_constraint_to_inner(md, **self._inner_constraint) # if _inner_constraint is not set, this does not change the md.
        return md
   
    def get_current_boxsize(self):
        """Get current box size in Angstrom."""
        boxsize = self.topology.getUnitCellDimensions()[0]._value  # assuming cubic box
        return boxsize 
    
    def warmup_simulation(self, time_steps=[0.5,1,2], steps=[10,50,1000], temperatures=[1,10,None], trajfilename_prefix='warmup', traj_frequency=10):
        """Perform stepwise warm-up MD.
        Args:
            time_steps (list): list of time steps in fs for each stage
            steps (list): list of number of steps for each stage
            temperatures (list): list of temperatures in K for each stage
            trajfilename_prefix (str): prefix for output DCD trajectory filenames, is combined with temperature of the stage
            traj_frequency (int): frequency of saving frames to DCD trajectory
        Returns:
            None
        """
        if temperatures[-1] is None:
            temperatures[-1] = self.default_temperature
            
        for dt, steps, temp in zip(time_steps, steps, temperatures):
            md = self.create_MD(
                dt=dt,
                equil_temp=temp,
                cv_atoms=self.indices_inner,
            )
            if self.save_traj:
                md.simulation.reporters.append(app.DCDReporter(f'{trajfilename_prefix}_{temp}K.dcd', traj_frequency))
            md.run(steps)
            self.positions = md.simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    
    def minimize_energy(self, tolerance=Quantity(10, unit.kilojoules_per_mole/unit.nanometer), max_iterations=1000):
        """Minimize energy of the system.
        Args:
            **kwargs: additional arguments for OpenMM minimizeEnergy() ()
        """
        md = self.create_MD(
            dt=self.default_dt,
            equil_temp=self.default_temperature,
            cv_atoms=self.indices_inner,
        )
        md.simulation.minimizeEnergy(tolerance=tolerance, maxIterations=max_iterations)
        self.positions = md.simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
        
    def save_positions(self, outfile):
        """Save current positions to PDB file.
        Args:
            outfile (str): output PDB filename
        """
        with open(outfile, 'w') as file:
            app.outfile.writeFile(self.topology, self.positions, file)
    
    def equilibrate_box(
        self,
        target_pressure=1.0* unit.atmosphere, 
        datafilename="nptsim.csv",
        numsteps_per_NPT=10000,
        logging_frequency=100,
        max_NPT_cycles=10,
        volume_threshold=1.3,
        density_threshold=0.005,
        trajfilename='trajBarostat.dcd',
        traj_frequency=10,
    ):
        """Conduct NPT equilibration using MonteCarloBarostat until box size and density converge.
        Args:
            target_pressure (Quantity): target pressure in atmosphere
            datafilename (str): output CSV filename for NPT data
            numsteps_per_NPT (int): number of steps per NPT cycle
            logging_frequency (int): frequency of logging in steps
            max_NPT_cycles (int): maximum number of NPT cycles
            volume_threshold (float): threshold for volume standard deviation for convergence in nm^3
            density_threshold (float): threshold for density standard deviation for convergence in g/mL
            trajfilename (str): output DCD filename for NPT trajectory
            traj_frequency (int): frequency of saving frames to DCD trajectory
        Returns:
            None    
            """
        mcb_force_index = self.system.addForce(MonteCarloBarostat(target_pressure, self.default_temperature * unit.kelvin, 25)) 
        md = self.create_MD(
            dt=self.default_dt,
            equil_temp=self.default_temperature,
            cv_atoms=self.indices_inner,
        )
        
        # To track positions 
        if self.save_traj:
            md.simulation.reporters.append(app.DCDReporter(trajfilename, traj_frequency))
 
        # create data file reporter 
        statedatareporter_file=app.StateDataReporter(datafilename, logging_frequency, step=True, time=True,
                                                        potentialEnergy=True, kineticEnergy=True, volume=True,
                                                        density=True, temperature=True, separator=',', append=True)
        md.simulation.reporters.append(statedatareporter_file)

        # create logging file
        with open(datafilename, 'w') as file:
            file.write('#"Step","Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"\n')

        numpoints_for_convergence_check=numsteps_per_NPT//logging_frequency

        # Starting parameters
        volume_std, density_std = 10, 1
        steps = 0

        # define function to read volume and density from the tracked data file
        def read_NPT_statefile(npt_output):
            import csv
            from collections import defaultdict
            # Read in CSV file of last NPT simulation and store in lists
            columns = defaultdict(list)

            with open(npt_output, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    for (k, v) in row.items():
                        columns[k].append(v)
            # Extract step number, volume and density and cast as floats
            steps = np.array(columns['#"Step"'])
            volume = np.array(columns["Box Volume (nm^3)"]).astype(float)
            density = np.array(columns["Density (g/mL)"]).astype(float)

            resultdict = {"steps": steps, "volume": volume, "density": density}
            return resultdict
        
        
        for i in range(max_NPT_cycles):
            md.run(numsteps_per_NPT)
            steps += numsteps_per_NPT
            
            NPTresults = read_NPT_statefile(datafilename)
            volume = NPTresults["volume"][-numpoints_for_convergence_check:]
            density = NPTresults["density"][-numpoints_for_convergence_check:]
            volume_std = np.std(volume)
            density_std = np.std(density)

            print("Total steps taken:", steps)
            print(f"Total simulation time: {md.dt/ 1000 * steps} ps")
            print("Current Volume:", volume[-1])
            print(f"Current Density: {density[-1]}")
            print()
            print(f"Current Volume SD: {volume_std}   (threshold: {volume_threshold})")
            print(f"Current Density SD: {density_std} (threshold: {density_threshold})")

            if volume_std < volume_threshold and density_std < density_threshold:
                print(f"Equilibration of periodic box finished after {steps} and {md.dt/ 1000 * steps} ps !\n")
                break

            if i == max_NPT_cycles-1:
                print(f"Warning: Max NPT cycles reached ({max_NPT_cycles}). Total steps taken: {steps} and {md.dt/ 1000 * steps} ps !\n")
                print("Warning: the NPT simulation may not be properly converged")
                break
        state = md.simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        self.positions = state.getPositions()
        # write out box vectors to topology to save new box size
        self.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        self.system.removeForce(mcb_force_index)
        
        
        