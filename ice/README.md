Proton-disordered ice structures
================================

This folder contains structures of proton-disordered hexagonal ice,
converted from the supporting information of Hayward and Reimers, JCP 1997.
Proton disorder indicates the random arrangement of hydrogen-bonding
(and hence molecular orientations) in a tetrahedral ice network.
Different orientations lead to subtle energy differences, and so the
energetics of different proton patterns is a very sensitive indicator
of rotational invariance. 

This folder contains uses i-PI with the PET driver, and contains two scripts:
`run_ice.sh` loops over the 9 orthorhombic ice supercells, and runs for 
each a geometry optimization. The `input-geop.xml` i-PI input file 
contains the parameters of the optimizer, and specifies the output files.
The script should be ran twice, once to run the optimization with 
a single PET model, once using rotational augmentation.

The `pproc.py` script analyzes the outputs, and extracts diagnostics of 
the differences between the "direct" and the "augmented" runs, such as 
differences in energies and forces, and difference in relaxed geometries.
