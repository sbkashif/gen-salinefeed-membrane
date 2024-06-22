# README

## Project Description

This project aims to ionize the feed side in a membrane-water systems to study salt trasnport through the membrane. The 'ionization' process entails modifying the `gro` and `top` files to add ions to the system. The ions are added to the feed side of the membrane only. The `gmx genion` tool available with gromacs does not support add ions in specific z-range (feed size in our case). Therefore, these customized scripts had to be written.


## Other details

The ionization is done such that the ion concentration is same as what was reported in Malmir and HajiAkbari, 2019. The ions are added to the feed side of the membrane only. The ions are added in customized z-range. The keep the system size miminal, the membrane is shifted downwards in the box so that miminal water can be added there while it is still enough to have the experimental water density.

## Logic 

- First, the gro file is trimmed to contain only the polymeric membrane and water molecules within the polymer, as determined from extremes of z-coordinates of the polymer. This is done using `obtain_polymer_domain`. By default, the script gives out the `cut.gro` file.
- Then, the `resolvate_and_salinate.sh` script has shell commands that operates on`cut.gro` to adjust the box such that it has target (feed ion concentration) z-dimensions of feedside and permeate side.
- Finally, `add_ions` (available in `domain_ionization` repo) and `modify_topology` are used to add ions to the feed side of the membrane and modify the topology file to include the ions. The ions are added in the feed side of the membrane only.

## Running the code
- Frist run the Jupyter notebook `determine_global_parameters.ipynb` to determine the global parameters for the system. The `PolymerDomain` class is called in the notebook, does all the caculations of finding feed side and permeate side of the membrane. The users are encourage to go through the table generated and ensure that the domain sizes are not drastically different across configurations.

- Run `resolvate_and_salinate.sh` to run the scripts dealing with modifying the `gro` and `itp` files. There is a manual entry for the variables `total_water_mol` and `ions_each_type` which is calculated in the notebook and has to be  entered in the shell script. In the future, this can be automated.