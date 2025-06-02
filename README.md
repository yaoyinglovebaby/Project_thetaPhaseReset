# Project_thetaPhaseReset
# Theta Phase Reset Analysis

This repository contains all code used for simulating and analyzing **theta phase reset** mechanisms in single-neuron spiking data. The goal of this project is to distinguish between phase reset and evoked response models of neural oscillations using simulations, spike-field coupling analysis, and statistical classification methods.

## Contents

- `simulation/`: MATLAB code for simulating spike trains and LFPs under:
  - Phase Reset (PR) model
  - Evoked Response (ERP) model
- `analysis/`: scripts for:
  - Oscillation score (O-score) calculation
  - Spike-LFP phase coupling (PPC) analysis
  - ROC curve generation and FDR correction
- `figures/`: example figures generated from simulations and analysis
- `utils/`: helper functions for signal processing and plotting

## Getting Started

1. Make sure you have MATLAB installed (tested with R2023a+).
2. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/theta-phase-reset.git
   cd theta-phase-reset
