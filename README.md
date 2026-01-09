# `createsys`: Create molecular systems for QM/MM hyperreactor dynamics
This package bundles some utilities for the fast setup of QM/MM hyperreactor simulations run with the [`adaptive_sampling`](https://github.com/ochsenfeld-lab/adaptive_sampling/tree/exploration_tools_qmmm) package. 
At some point, these packages might be merged.

## Available features
* **intial conformation generation** of a mixture of reactants in a sphere sourrounded by a cubical box of solvent molecules (by accessing [`packmol`](https://m3g.github.io/packmol/) through a simple wrapper). For each molecule, a single `.pdb` is required.

* **equilibration routines** based on `openmm`/`adaptive_sampling`:

    * energy minimization

    * multi-step warmup MD

    * Box equilibration MD with Monte-Carlo-Barostat

> [!NOTE]
> A spherical confinement can be added to equilibration MDs to keep the reactants centered.

* automated **force field setup** based on a single `.pdb`. With geometry-inferred bond orders, an `openff`-SMIRNOFF-based force field (`openff-2.2.1.offxml`) is created together with a suitable topology. This function prioritizes convenience over customizability.

## Examples
For initial configuration generation of a non-enzymatic nucleoside synthesis example, see [here](examples/initial_configuration.py). 

Equilibration routines are found [here](examples/equilbration.py).

For the topology and force field generation, simply use
```python
from createsys import createTopologyAndForceFieldFromPDB

topology, forcefield = createTopologyAndForceFieldFromPDB("mystructure.pdb")
```
The connectivity information of the PDB is used to determine molecular connectivity. Please note that bond orders are inferred from the geometry.

To identify the inner region of a structure created with `createsys`, you can use
```python
from createsys import getAtomIndicesFromResidueNames
atom_indices_inner_sphere = getAtomIndicesFromResidueNames("mystructure.pdb", ["INN"])

```
## Installation
As this package is an add-on to the `adaptive_sampling` package, please follow the instructions given [here](https://github.com/ochsenfeld-lab/adaptive_sampling/blob/exploration_tools_qmmm/README.md#install). This also installs `createsys`.