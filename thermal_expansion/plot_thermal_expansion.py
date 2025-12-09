#!/usr/bin/env python3
"""
Plot QHA thermal expansion results from megafit.
Output file has hexagonal a and c already in Angstrom.
"""

import numpy as np
import matplotlib.pyplot as plt

# Read the equilibrium data
data = np.loadtxt('/Users/ds.kim/Software/tdep/thermal_expansion/outfile.ac_equilibrium.dat')
T = data[:, 0]       # Temperature (K)
a_eq = data[:, 1]    # a lattice parameter (Angstrom) - hexagonal
c_eq = data[:, 2]    # c lattice parameter (Angstrom) - hexagonal
F_min = data[:, 3]   # Free energy minimum (eV/atom)

# These are already in Angstrom (hexagonal cell parameters)
a_angstrom = a_eq
c_angstrom = c_eq

# Reference values at lowest temperature
a0 = a_eq[0]
c0 = c_eq[0]

# Calculate thermal expansion (relative change)
delta_a = (a_eq - a0) / a0 * 100  # percent
delta_c = (c_eq - c0) / c0 * 100  # percent

# Volume (for hexagonal: V = (sqrt(3)/2) * a^2 * c)
V = (np.sqrt(3) / 2) * a_eq**2 * c_eq
V0 = V[0]
delta_V = (V - V0) / V0 * 100  # percent

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. Free energy vs temperature
ax1 = axes[0, 0]
ax1.plot(T, F_min, 'b-', linewidth=2)
ax1.set_xlabel('Temperature (K)', fontsize=12)
ax1.set_ylabel('Free Energy (eV/atom)', fontsize=12)
ax1.set_title('QHA Free Energy vs Temperature', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 1000)

# 2. Lattice parameters vs temperature
ax2 = axes[0, 1]
ax2.plot(T, a_angstrom, 'b-', linewidth=2, label='a')
ax2.set_xlabel('Temperature (K)', fontsize=12)
ax2.set_ylabel('a (Å)', fontsize=12, color='b')
ax2.tick_params(axis='y', labelcolor='b')
ax2.set_xlim(0, 1000)

ax2b = ax2.twinx()
ax2b.plot(T, c_angstrom, 'r-', linewidth=2, label='c')
ax2b.set_ylabel('c (Å)', fontsize=12, color='r')
ax2b.tick_params(axis='y', labelcolor='r')

ax2.set_title('Lattice Parameters vs Temperature', fontsize=14)
ax2.grid(True, alpha=0.3)

# Add legend
lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2b.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# 3. Thermal expansion (relative change)
ax3 = axes[1, 0]
ax3.plot(T, delta_a, 'b-', linewidth=2, label='Δa/a₀')
ax3.plot(T, delta_c, 'r-', linewidth=2, label='Δc/c₀')
ax3.plot(T, delta_V, 'k--', linewidth=2, label='ΔV/V₀')
ax3.set_xlabel('Temperature (K)', fontsize=12)
ax3.set_ylabel('Relative Change (%)', fontsize=12)
ax3.set_title('Thermal Expansion', fontsize=14)
ax3.legend(loc='upper left')
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 1000)
ax3.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)

# 4. Thermal expansion coefficients (derivative)
ax4 = axes[1, 1]
# Calculate thermal expansion coefficient: alpha = (1/L) * dL/dT
# Use finite differences
dT = np.diff(T)
alpha_a = np.diff(a_eq) / a_eq[:-1] / dT * 1e6  # in 10^-6 K^-1
alpha_c = np.diff(c_eq) / c_eq[:-1] / dT * 1e6
alpha_V = np.diff(V) / V[:-1] / dT * 1e6
T_mid = (T[:-1] + T[1:]) / 2

ax4.plot(T_mid, alpha_a, 'b-', linewidth=2, label='αₐ')
ax4.plot(T_mid, alpha_c, 'r-', linewidth=2, label='αc')
ax4.plot(T_mid, alpha_V / 3, 'k--', linewidth=2, label='αᵥ/3')
ax4.set_xlabel('Temperature (K)', fontsize=12)
ax4.set_ylabel('Thermal Expansion Coeff. (10⁻⁶ K⁻¹)', fontsize=12)
ax4.set_title('Thermal Expansion Coefficients', fontsize=14)
ax4.legend(loc='upper right')
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 1000)
ax4.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)

plt.tight_layout()
plt.savefig('/Users/ds.kim/Software/tdep/thermal_expansion/qha_thermal_expansion.png', dpi=150, bbox_inches='tight')
plt.show()

print("Saved: qha_thermal_expansion.png")

# Print some key values
print("\n" + "="*60)
print("Summary of QHA Results for LiNbO₃ (hexagonal cell)")
print("="*60)
print(f"\nAt T = 0 K:")
print(f"  a = {a_eq[0]:.4f} Å")
print(f"  c = {c_eq[0]:.4f} Å")
print(f"  c/a = {c_eq[0]/a_eq[0]:.4f}")
print(f"  F = {F_min[0]:.6f} eV/atom")

print(f"\nAt T = 300 K (interpolated):")
idx_300 = np.argmin(np.abs(T - 300))
print(f"  a = {a_eq[idx_300]:.4f} Å")
print(f"  c = {c_eq[idx_300]:.4f} Å")
print(f"  c/a = {c_eq[idx_300]/a_eq[idx_300]:.4f}")
print(f"  Δa/a₀ = {delta_a[idx_300]:.4f} %")
print(f"  Δc/c₀ = {delta_c[idx_300]:.4f} %")

print(f"\nAt T = 1000 K:")
print(f"  a = {a_eq[-1]:.4f} Å")
print(f"  c = {c_eq[-1]:.4f} Å")
print(f"  c/a = {c_eq[-1]/a_eq[-1]:.4f}")
print(f"  Δa/a₀ = {delta_a[-1]:.4f} %")
print(f"  Δc/c₀ = {delta_c[-1]:.4f} %")
