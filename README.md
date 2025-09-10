# hertault-model (HCM) 
# 🌌 Hertault Model: Dark Sector Unification at the Cosmological Critical Density

![Hertault Model Results](figures/hertault_main_results.png)

Welcome to the official repository for the paper **"Dark Sector Unification at the Cosmological Critical Density"** by Hugo Hertault (Lille, France, September 2025). This repository contains the complete numerical simulation suite that produced the results presented in the paper, proposing a novel unified theory where dark matter and dark energy emerge from a single scalar field φ with environment-dependent dynamics.

## 📝 Overview

The Hertault model introduces a scalar field φ that undergoes a natural phase transition at the cosmological critical density \( \rho_c \approx 3.64 \times 10^{-47} \, \text{GeV}^4 \), behaving as:
- **Dark Matter** (\( \langle w_\phi \rangle \simeq 0 \)) at high density (\( \rho_m \gg \rho_c \)).
- **Dark Energy** (\( w_\phi \simeq -1 \)) at low density (\( \rho_m \ll \rho_c \)).

### Key Achievements
- **Resolves the cosmic coincidence problem** without fine-tuning, with a transition at \( z_{\text{trans}} \approx 0.45 \).
- **Reduces cosmological tensions**:
  - **H₀ tension**: 60% reduction (from 5.4σ to 2.1σ).
  - **σ₈ tension**: 50% reduction (from 2.6σ to 1.3σ).
- **Addresses the neutrino mass puzzle** through environment-dependent neutrino masses.
- **Satisfies all experimental constraints** with a universal coupling \( \alpha = 7.8 \times 10^{-4} \).
- **Predicts testable signatures** for DESI, Euclid, Einstein Telescope, and MICROSCOPE-2.

The code is designed to be accessible, modular, and reproducible, with publication-quality figures and validation tests.

## 🚀 Features

- **HertaultField**: Implements the scalar field with an environment-dependent potential and universal coupling.
- **CosmologicalEvolution**: Solves the coupled field-Friedmann equations for background evolution.
- **PerturbationEvolution**: Computes linear perturbations and modified growth factor.
- **ObservationalConstraints**: Analyzes cosmological tensions (H₀, σ₈) against Planck, SH0ES, KiDS, and DES data.
- **ParameterOptimization**: Optimizes model parameters using differential evolution.
- **Publication Figures**: Generates high-quality plots, including phase diagrams, field evolution, Hubble parameter evolution, and tension reduction.

## 📊 Key Results

- **H₀ Prediction**: \( 69.8 \pm 1.2 \, \text{km/s/Mpc} \)
- **σ₈ Prediction**: \( 0.789 \pm 0.011 \)
- **Transition Redshift**: \( z_{\text{trans}} \approx 0.45 \pm 0.05 \)
- **Equation of State**: \( w_\phi(z=0) \approx -0.97 \pm 0.02 \)
- **Tension Reductions**:
  - H₀: Reduced from 5.4σ to 2.1σ (61% improvement).
  - σ₈: Reduced from 2.6σ to 1.3σ (50% improvement).
- **Experimental Constraints**:
  - Fifth force: \( \alpha < 10^{-3} \) ✓
  - Equivalence principle (MICROSCOPE): \( \eta < 10^{-13} \) ✓
  - Gravitational wave speed: \( |c_{\text{gw}} - c|/c < 10^{-15} \) ✓
  - Stellar cooling: \( \alpha < 5 \times 10^{-3} \) ✓

## 🛠 Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/hugohertault/hertault-model.git
   cd hertault-model
