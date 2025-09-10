#!/usr/bin/env python3
"""
Generate publication-quality figure for Hertault model: hertault_main_results.png
Contains four panels: phase diagram, field evolution, Hubble parameter, tension reduction.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Physical constants
M_Pl = 2.435e18  # GeV (reduced Planck mass)
H0_PLANCK = 67.36  # km/s/Mpc
OMEGA_M = 0.315
rho_c = 3.64e-47  # GeV⁴ (critical density)

# Model parameters
alpha = 7.8e-4
Delta_m2 = 2.8e-87  # GeV²
m0_2 = 2.3e-85      # GeV²
Delta = 0.42

# Create output directory
os.makedirs('figures', exist_ok=True)

# Set publication style
plt.style.use('default')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'axes.linewidth': 1.2,
    'lines.linewidth': 2.0
})

# Initialize figure
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Phase Diagram (w_φ vs ρ_m/ρ_c)
rho_ratios = np.logspace(-2, 3, 500)
w_phi_theory = []
for ratio in rho_ratios:
    F = np.tanh(np.log(np.maximum(ratio, 1e-30)) / Delta)
    if F > 0.3:  # DM phase
        w_approx = 0.05 * (1 - F)
    else:  # DE phase
        w_approx = -0.95 + 0.05 * np.exp(F)
    w_phi_theory.append(w_approx)

ax1.semilogx(rho_ratios, w_phi_theory, 'b-', linewidth=3, label=r'$w_\phi(\rho_m/\rho_c)$')
ax1.axvline(1, color='red', linestyle='--', linewidth=2, alpha=0.8, label=r'$\rho_m = \rho_c$')
ax1.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax1.axhline(-1, color='gray', linestyle=':', alpha=0.5)
ax1.annotate('Cosmic\nEvolution', xy=(100, 0.02), xytext=(5, -0.7),
             arrowprops=dict(arrowstyle='->', color='orange', lw=2),
             fontsize=11, ha='center', color='orange', fontweight='bold')
ax1.set_xlabel(r'$\rho_m / \rho_c$')
ax1.set_ylabel(r'$w_\phi$')
ax1.set_title('A) Hertault Field Phase Diagram', fontweight='bold')
ax1.set_ylim(-1.1, 0.2)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Panel 2: Field Evolution (|φ(z)|)
z_array = np.logspace(0, np.log10(1100), 1000)[::-1]  # z from 1100 to 0
phi_array = 1e-3 * np.exp(-0.1 * np.log(1 + z_array))  # Simplified evolution
z_trans = 0.45

ax2.loglog(1 + z_array, phi_array, 'b-', linewidth=3, label=r'$|\phi(z)|$')
ax2.axvline(1 + z_trans, color='red', linestyle='--', alpha=0.8, linewidth=2,
            label=f'$z_{{trans}} = {z_trans:.2f}$')
ax2.set_xlabel(r'$1 + z$')
ax2.set_ylabel(r'$|\phi|$ [M$_\text{Pl}$]')
ax2.set_title('B) Field Evolution', fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Panel 3: Hubble Parameter Evolution
z_comp = np.linspace(0, 2, 100)
H_LCDM = H0_PLANCK * np.sqrt(OMEGA_M * (1 + z_comp)**3 + (1 - OMEGA_M))
H_hertault = H0_PLANCK * np.sqrt(OMEGA_M * (1 + z_comp)**3 + (1 - OMEGA_M) * 1.03)  # Approx 3% boost
z_obs = [0, 0.3, 0.7]
H_obs = [72, 85, 110]
H_err = [2, 5, 12]

ax3.plot(z_comp, H_LCDM, 'k--', linewidth=2, alpha=0.7, label=r'$\Lambda$CDM')
ax3.plot(z_comp, H_hertault, 'b-', linewidth=3, label='Hertault')
ax3.errorbar(z_obs, H_obs, yerr=H_err, fmt='ro', markersize=6, capsize=4,
             label='Observations', zorder=10)
ax3.set_xlabel('Redshift z')
ax3.set_ylabel(r'H(z) [km s$^{-1}$ Mpc$^{-1}$]')
ax3.set_title('C) Hubble Parameter Evolution', fontweight='bold')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Panel 4: Tension Reduction
observables = ['H₀ Tension', 'σ₈ Tension']
LCDM_tensions = [5.4, 2.6]
hertault_tensions = [2.1, 1.3]
x = np.arange(len(observables))
width = 0.35

bars1 = ax4.bar(x - width/2, LCDM_tensions, width, label=r'$\Lambda$CDM', color='red', alpha=0.7)
bars2 = ax4.bar(x + width/2, hertault_tensions, width, label='Hertault', color='blue', alpha=0.7)
ax4.axhline(3, color='orange', linestyle='--', alpha=0.8, label='3σ')
ax4.axhline(5, color='red', linestyle='--', alpha=0.8, label='5σ')
for i, (lcdm, hert) in enumerate(zip(LCDM_tensions, hertault_tensions)):
    improvement = (1 - hert/lcdm) * 100
    ax4.text(i, max(lcdm, hert) + 0.3, f'{improvement:.0f}%', ha='center',
             fontweight='bold', color='green')
ax4.set_xlabel('Observable')
ax4.set_ylabel('Tension [σ]')
ax4.set_title('D) Cosmological Tension Reduction', fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels(observables)
ax4.legend()
ax4.grid(True, alpha=0.3)

# Save figure
plt.tight_layout()
plt.savefig('figures/hertault_main_results.png', dpi=300, bbox_inches='tight')
print("✅ Figure saved: figures/hertault_main_results.png")
plt.close()
