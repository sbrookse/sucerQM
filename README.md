# SUCER-qm 
### Saccharide pUcker Conformation and Energy exploreR using quantum mechanics

![Example for xylose](xylose.gif)

## Overview
- sucer-QM includes two Python packages (pucker_gen and pucker_analysis) that enable the generation and energetic exploration of puckered conformations for any 5-membered (furanose) or 6-membered (pyranose) substituted ring. While the program is designed with monosaccharides in mind, it is compatible with almost all 5- and 6- membered rings.
- Accepts SMILES strings as input and generates puckered conformation files in pdb/xyz or Gausian input format based on specified Cremer-Pople puckering parameters (q, $\theta$ for furanoses and q, $\theta$, $\phi$ for pyranoses)
- Ensuring energy minima for each pucker state involves the full exploration of the orientation of all exocyclic groups attached to the monosaccharide. To enable this, sucer-QM has pre-defined functions (starting with inital_dihedral_scan) that enable QM energetic scans of dihedral angles of exocyclic groups. This is undertaken in an iterative manner with (repeat_dihedral_scan) functions, wherein each exocylcic group is scanned independently until no lower energy minima are observed.
- Finally, the final_unconstrained  or final_constrained functions enable the calculation of final optimized structures for thermodynamic and energetic properties.

 Contributor: [Sean Brooks](https://github.com/sbrookse)

## Dependencies

This code requires the following pre-requisite packages to be installed.
- RDKit
- numpy
- pandas
- Matplotlib
- RingAnalysis (https://github.com/lucianlschan/RING)
- RingReconstruction (https://github.com/lucianlschan/RING)
- py_rdl (RingDecomposerLib) (https://github.com/rareylab/RingDecomposerLib)


## Running sucer-QM
* Step 0
  - Update the SMILES strings of the sugar molecules whose pucker states are desired to be explored in the Sugar_Info.xlsx spreadsheet. Type "yes" in the Include field for a given sugar to include that sugar in the current run and type "no" if it should not be included.
  - Edit the c5_pucker_params.csv and/or c6_pucker_params.csv file(s) to ensure that the parameters (esp. q values) correspond to the desired pucker states to be explored. Previously validated parameters from Chan, et al., (2020) are used as default in the spreadsheet but can be changed depending on user needs.
* Step 1 - Run initial optimizations for each sugar and its various pucker states
  - ```
    python
    import pucker_gen as pg
    pg.generate_puckers('input')
    ```
    - This will generate a folder called input containing folders named after each sugar whose pucker state is desired to be explored. If using relative path, the folder will generate in the same folder where the pucker_gen.py script is saved; this can be changed as desired using an absolute file path.
    - Each sugar folder contains sub-folders named after the pucker state and contains Gaussian 16 input files (.com) and .pdb files for the pucker state.
  - If you have access to an HPC workstation with Gaussian 16 installed, then you may run the following command to submit the QM calcualtions
  - ```
    sh submit.sh
    sh compile.sh
    ```
    - 
* Step 2 - Perform dihedral scans for exocyclic groups attached to the sugar
  - kdjfa'lsdf
  - ```
    ksjdfasdf
    asdfasdf
    ```
    - This will generate gaussian input files for running dihedral scans for each individual exocyclic group.
    -
* Step 3 - Perform final optimizations of the lowest energy configurations identified from Step 2
  - asdfkjasdf
  - asdfkljasldkfj
  - 
aaaa
