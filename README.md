# Eco-Evolutionary Dynamics in Age-Structured Populations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
*(Note: Update the Zenodo badge with your actual DOI once generated)*

## Overview
This repository contains the source code, simulation scripts, and mathematical derivations accompanying the manuscript: **"Eco-Evolutionary Dynamics in Age-Structured Populations: A Hybrid Approach of Leslie-Lefkovitch Matrices and Density-Dependent Gamete Pools"** (currently under review).

Our framework resolves the mathematical tension between structured population ecology and Mendelian genetics. By decoupling the deterministic process of survival (via Leslie-Lefkovitch matrices) from the stochastic mechanics of sexual reproduction (via a simulated gamete pool), the model maintains a low-dimensional architecture ($s \times s$) while explicitly tracking allele frequencies and non-linear density-dependent regulation (Ricker).

## Key Features
* **Hybrid Architecture:** Integrates structured demography with single-locus Mendelian segregation without inflating the matrix dimensionality to a "megamatrix".
* **Density-Dependent Scramble Competition:** Utilizes a Ricker formulation applied prior to demographic transitions to simulate severe resource competition and demographic overshoots.
* **Transient & Stability Analysis:** Includes scripts for 1D sensitivity sweeps, 2D stability heatmaps (Coefficient of Variation), and local stability evaluations (Jacobian eigenvalues) to map the boundaries between stable equilibria, limit cycles, and deterministic chaos.

## Repository Structure
The repository is organized to ensure full reproducibility of the figures and analyses presented in the manuscript:

```text
├── src/                      # Core model functions and classes
│   ├── hybrid_model.py       # Main simulation algorithm (Gamete pool + Leslie)
│   └── stability_utils.py    # Functions for CV and Jacobian calculations
├── scripts/                  # Scripts to reproduce manuscript figures
│   ├── fig1_buteo.py         # Empirical validation using Buteo buteo vital rates
│   ├── fig2_archetype.py     # Fast-lived archetype transient dynamics
│   ├── fig3_1D_sweep.py      # One-dimensional survival advantage sweep
│   └── fig4_2D_heatmap.py    # 2D stability landscape and contour mapping
├── data/                     # (Optional) Any static parameter files or empirical vital rates
└── README.md


## Requirements
The code is written in Python 3.8+ and relies on standard scientific libraries. You can install the dependencies using:Bashpip install numpy scipy matplotlib

## Usage
To replicate the figures from the manuscript, simply run the corresponding script from the terminal. For example, to generate the 2D Stability Landscape (Figure 4):

Bash python scripts/fig4_2D_heatmap.py

Note: The 2D heatmap script simulates thousands of multigenerational evolutionary trajectories and may take a few minutes to complete depending on your hardware.

## Mathematical Appendix & Code

The script stability_utils.py contains the computational implementation of the analytical derivation described in Appendix B of the manuscript. It constructs the Jacobian matrix at the non-trivial equilibrium point ($$N^*$$) and calculates the dominant eigenvalue ($$|\lambda_{dom}|$$) to explicitly detect the period-2 bifurcation thresholds.

## Citation
If you use this code or framework in your research, please cite the original paper:Alcántara-Rodríguez, J. A. (2026). Eco-Evolutionary Dynamics in Age-Structured Populations: A Hybrid Approach of Leslie-Lefkovitch Matrices and Density-Dependent Gamete Pools.

## License
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details. 
You are free to use, modify, and distribute this software, but any derivative works must also be distributed under the same open-source terms.
