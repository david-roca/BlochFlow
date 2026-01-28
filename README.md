# BlochFlow

**A MATLAB solver for Bloch analysis of flow perturbations in a channel over a lattice of phononic subsurface units**

This repository contains a spectral-FEM solver that performes Bloch analysis of flow perturbations in a channel with a lattice of phononic subsurface units. The code computes the Bloch dispersion relations (frequency vs. wavenumber) and mode shapes characterizing the stability characteristics of a channel flow interacting with a periodic arrangement of phononic subsurface (PSub) units. It considers a single unit cell and solves the generalized Orr-Sommerfeld equation coupled with an elastodynamics model for the surface admittance representing the PSub motion.

![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-orange)

## Features

* **Hybrid discretization:** Uses Fourier modes in the streamwise and spanwise directions and Finite Elements (FE) in the wall-normal direction.
* **Fluid-Structure Interaction:** Fully coupled formulation linking the Orr-Sommerfeld fluid equation with a 1D PSub solid model via impedance boundary conditions.
* **Metamaterial support:** Capable of modeling standard Phononic Crystals (PnC) and Locally Resonant Metamaterials (MM) with internal resonators.
* **Flexible geometry:** Supports various lattice shapes (square, hexagon).

## Installation

1.  Clone the repository:
```bash
git clone [https://github.com/david-roca/BlochFlow.git](https://github.com/david-roca/BlochFlow.git)
```
2.  Open MATLAB and navigate to the folder.
3.  Add the source code to your path:
```matlab
addpath(genpath('src'));
savepath;
```

## Usage

A complete example is provided in `main_example.m`.

## Citation

If you use this solver in your work, please cite the following paper:
> **M.I. Hussein et al. (2025).** "Scatterless interferences: Delay of laminar-to-turbulent flow transition by a lattice of subsurface phonons". *Preprint*. DOI: 10.48550/arXiv.2503.18835
```bibtex
@article{Hussein2025,
  author = {Hussein, Mahmoud I. and Roca, David and Harris, Adam R. and Kianfar, Armin},
  title = {Scatterless interferences: Delay of laminar-to-turbulent flow transition by a lattice of subsurface phonons},
  journal = {Preprint},
  year = {2025},
  doi = {10.48550/arXiv.2503.18835}
}
```
