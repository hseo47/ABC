# ABC-Adhesion-of-Bacteria-to-Chitosan

This repository contains MATLAB code and Scanned Electron Microscopy (SEM) image files used in the study:  
**Enhancing Bacterial Adhesion with Hydro-Softened Chitosan Film**  
The code simulates bacterial adhesion using a kinetic Monte Carlo framework

## Contents

- 'kMCsimulation.m' – Runs the stochastic kMC simulation loop (JKR-Griffith adhesion threshold, xDLVO)
- 'morphImgPro.m' – Runs the morphological classification of SEM images post-incubation
- 'Data.rtf' – Morphological classification outputs
- 'Simulatd Data.rtf' - Simulation outputs (Bacteria Count) over 24 hr run
- Simulated binary grid (kMC) and SEM images

## Requirements

- MATLAB R2023 or later (Preferably)
- Image Processing Toolbox

Tested with MATLAB_R2024b; No additional external libraries are needed

## Usage 

- Run MATLAB
- Run kMCsimulation.m (Adjust parameters if necessary)
- Run morImgPro.m (Change input file name)

## Code Availability Statement 

In compliance with institutional data management policies, the full codebase is maintained on Georgia Tech’s GitHub Enterprise instance. A stable public version, sufficient for reproducing the results reported in this manuscript, is available at https://github.com/hseo47/ABC_Adhesion_of_Bacteria_to_Chitosan and has been archived with a permanent DOI at https://doi.org/10.5281/zenodo.15323528. This public version includes usage instructions, example datasets, and all necessary dependencies for replication.
