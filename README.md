# bacterial-adhesion-kMC

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

The repository will be made public upon publication
