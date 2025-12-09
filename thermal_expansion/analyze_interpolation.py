#!/usr/bin/env python3
"""
Analyze megafit interpolation scheme and plot free energies at training points.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

print("="*70)
print("MEGAFIT INTERPOLATION SCHEME ANALYSIS")
print("="*70)

#------------------------------------------------------------------------------
# 1. Parse infile.simulations
#------------------------------------------------------------------------------
print("\n1. INFILE.SIMULATIONS FORMAT:")
print("-"*70)

with open('/Users/ds.kim/Software/tdep/thermal_expansion/infile.simulations', 'r') as f:
    lines = f.readlines()

ndim = int(lines[0].strip())
nsim = int(lines[1].strip())
dim_names = lines[2].strip().split()
order_per_dim = [int(x) for x in lines[3].strip().split()]
reference = lines[4].strip()

print(f"  Number of dimensions:     {ndim}")
print(f"  Number of simulations:    {nsim}")
print(f"  Dimension names:          {dim_names}")
print(f"  Polynomial order per dim: {order_per_dim}")
print(f"  Reference structure:      {reference}")

# Parse simulation data
a_vals = []
c_vals = []
energies = []
for i in range(5, 5 + nsim):
    parts = lines[i].split()
    a_vals.append(float(parts[0]))
    c_vals.append(float(parts[1]))
    energies.append(float(parts[2]))

a_vals = np.array(a_vals)
c_vals = np.array(c_vals)
energies = np.array(energies)

# Convert Bohr to Angstrom for display
bohr_to_angstrom = 0.529177249

print(f"\n  Training grid:")
print(f"    a range: {a_vals.min():.4f} to {a_vals.max():.4f} Bohr")
print(f"           = {a_vals.min()*bohr_to_angstrom:.4f} to {a_vals.max()*bohr_to_angstrom:.4f} Å")
print(f"    c range: {c_vals.min():.4f} to {c_vals.max():.4f} Bohr")
print(f"           = {c_vals.min()*bohr_to_angstrom:.4f} to {c_vals.max()*bohr_to_angstrom:.4f} Å")
print(f"    Static energy range: {energies.min():.6f} to {energies.max():.6f} eV/atom")

#------------------------------------------------------------------------------
# 2. Parse infile.evalpoints
#------------------------------------------------------------------------------
print("\n\n2. INFILE.EVALPOINTS FORMAT:")
print("-"*70)

with open('/Users/ds.kim/Software/tdep/thermal_expansion/infile.evalpoints', 'r') as f:
    eval_lines = f.readlines()

evalmode = int(eval_lines[0].strip())
eval_modes = {
    1: "Volume grid (1D)",
    2: "Point list",
    3: "Volume-Temperature grid (2D)",
    4: "a-c lattice parameter grid (2D)",
    5: "a-c-T grid (3D)"
}

print(f"  Eval mode: {evalmode} = {eval_modes.get(evalmode, 'Unknown')}")

if evalmode == 4:
    a_grid_info = eval_lines[1].strip().split()
    a_range = eval_lines[2].strip().split()
    c_grid_info = eval_lines[3].strip().split()
    c_range = eval_lines[4].strip().split()
    
    print(f"  a-axis grid: {a_grid_info[0]} points, {a_grid_info[1]} spacing")
    print(f"    Range: {a_range[0]} to {a_range[1]} Bohr")
    print(f"  c-axis grid: {c_grid_info[0]} points, {c_grid_info[1]} spacing")
    print(f"    Range: {c_range[0]} to {c_range[1]} Bohr")

#------------------------------------------------------------------------------
# 3. Explain interpolation orders
#------------------------------------------------------------------------------
print("\n\n3. INTERPOLATION ORDERS USED IN MEGAFIT:")
print("-"*70)

print("""
  ┌─────────────────────────────────────────────────────────────────────┐
  │                    MULTI-LEVEL INTERPOLATION                        │
  ├─────────────────────────────────────────────────────────────────────┤
  │                                                                     │
  │  LEVEL 1: Force Constants (from infile.simulations)                 │
  │  ─────────────────────────────────────────────────────────────────  │
  │  • Order per dimension: from line 4 of infile.simulations           │
  │  • Your setting: order = {} for 'a', order = {} for 'c'             │
  │  • This means Φ_αβ(a,c) = Σ c_ij * a^i * c^j   (i,j ≤ 2)            │
  │  • Total FC polynomial terms: (2+1) × (2+1) = 9 terms               │
  │                                                                     │
  │  LEVEL 2: Static Energy U(a,c)                                      │
  │  ─────────────────────────────────────────────────────────────────  │
  │  • Fitted from DFT energies in infile.simulations (3rd column)      │
  │  • Uses 4th order 2D polynomial: a^i × c^j where i,j ∈ [0,4]        │
  │  • Total terms: 5 × 5 = 25 coefficients                             │
  │  • This is the cold (T=0) DFT energy surface                        │
  │                                                                     │
  │  LEVEL 3: Phonon Free Energy F_ph(a,c,T)                            │
  │  ─────────────────────────────────────────────────────────────────  │
  │  • Computed at each (a,c) point using force constants from Level 1  │
  │  • At each T: fit F_ph(a,c) with 4th order 2D polynomial            │
  │  • Total terms: 25 coefficients per temperature                     │
  │                                                                     │
  │  LEVEL 4: Total Free Energy F_tot(a,c,T)                            │
  │  ─────────────────────────────────────────────────────────────────  │
  │  • F_tot = U(a,c) + U0(a,c) + F_ph(a,c,T)                            │
  │  • Combined 4th order polynomial fit at each T                      │
  │  • 25 polynomial coefficients stored per temperature                │
  │                                                                     │
  │  MINIMUM FINDING                                                    │
  │  ─────────────────────────────────────────────────────────────────  │
  │  • At each T: find (a*,c*) that minimizes F_tot(a,c,T)               │
  │  • Uses Newton-Raphson on the polynomial surface                    │
  │  • Output: a*(T), c*(T), F_min(T) in outfile.ac_equilibrium.dat     │
  │                                                                     │
  └─────────────────────────────────────────────────────────────────────┘
""".format(order_per_dim[0], order_per_dim[1]))

#------------------------------------------------------------------------------
# 4. Plot static energies at training points
#------------------------------------------------------------------------------
print("\n4. PLOTTING FREE ENERGIES AT TRAINING POINTS...")
print("-"*70)

# Get unique a and c values
a_unique = np.sort(np.unique(a_vals))
c_unique = np.sort(np.unique(c_vals))
na = len(a_unique)
nc = len(c_unique)

# Reshape energies into 2D grid
E_grid = np.zeros((na, nc))
for i, a in enumerate(a_unique):
    for j, c in enumerate(c_unique):
        idx = np.where((np.abs(a_vals - a) < 0.001) & (np.abs(c_vals - c) < 0.001))[0]
        if len(idx) > 0:
            E_grid[i, j] = energies[idx[0]]

# Create figure
fig = plt.figure(figsize=(16, 12))

# 1. 3D surface plot
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
A, C = np.meshgrid(a_unique, c_unique, indexing='ij')
surf = ax1.plot_surface(A * bohr_to_angstrom, C * bohr_to_angstrom, E_grid, 
                        cmap=cm.viridis, alpha=0.8)
ax1.scatter(a_vals * bohr_to_angstrom, c_vals * bohr_to_angstrom, energies, 
           c='red', s=50, label='Training points')
ax1.set_xlabel('a (Å)')
ax1.set_ylabel('c (Å)')
ax1.set_zlabel('E (eV/atom)')
ax1.set_title('Static DFT Energy Surface U(a,c)')
ax1.view_init(elev=20, azim=45)

# 2. Contour plot of static energy
ax2 = fig.add_subplot(2, 2, 2)
contour = ax2.contourf(A * bohr_to_angstrom, C * bohr_to_angstrom, E_grid, 
                       levels=30, cmap=cm.viridis)
ax2.scatter(a_vals * bohr_to_angstrom, c_vals * bohr_to_angstrom, 
           c='red', s=50, edgecolors='white', label='Training points')
ax2.set_xlabel('a (Å)')
ax2.set_ylabel('c (Å)')
ax2.set_title('Static Energy Contour U(a,c) [eV/atom]')
plt.colorbar(contour, ax=ax2, label='E (eV/atom)')

# Mark the minimum
min_idx = np.argmin(energies)
ax2.scatter(a_vals[min_idx] * bohr_to_angstrom, c_vals[min_idx] * bohr_to_angstrom, 
           c='white', s=200, marker='*', edgecolors='black', linewidths=2,
           label=f'Minimum: ({a_vals[min_idx]*bohr_to_angstrom:.4f}, {c_vals[min_idx]*bohr_to_angstrom:.4f}) Å')
ax2.legend()

# 3. Energy vs a (at different c values)
ax3 = fig.add_subplot(2, 2, 3)
colors = cm.plasma(np.linspace(0, 1, nc))
for j, c in enumerate(c_unique):
    idx = [i for i in range(len(c_vals)) if np.abs(c_vals[i] - c) < 0.001]
    ax3.plot(a_vals[idx] * bohr_to_angstrom, energies[idx], 'o-', 
            color=colors[j], label=f'c = {c*bohr_to_angstrom:.3f} Å')
ax3.set_xlabel('a (Å)')
ax3.set_ylabel('E (eV/atom)')
ax3.set_title('Energy vs a at fixed c')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)

# 4. Energy vs c (at different a values)
ax4 = fig.add_subplot(2, 2, 4)
colors = cm.plasma(np.linspace(0, 1, na))
for i, a in enumerate(a_unique):
    idx = [j for j in range(len(a_vals)) if np.abs(a_vals[j] - a) < 0.001]
    ax4.plot(c_vals[idx] * bohr_to_angstrom, energies[idx], 'o-', 
            color=colors[i], label=f'a = {a*bohr_to_angstrom:.3f} Å')
ax4.set_xlabel('c (Å)')
ax4.set_ylabel('E (eV/atom)')
ax4.set_title('Energy vs c at fixed a')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/ds.kim/Software/tdep/thermal_expansion/training_energies.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nSaved: training_energies.png")

#------------------------------------------------------------------------------
# 5. Summary table
#------------------------------------------------------------------------------
print("\n\n5. TRAINING POINT ENERGIES:")
print("-"*70)
print(f"{'a (Å)':>10} {'c (Å)':>10} {'E (eV/atom)':>15} {'ΔE (meV/atom)':>15}")
print("-"*70)

E_min = energies.min()
for i in range(len(energies)):
    dE = (energies[i] - E_min) * 1000  # convert to meV
    print(f"{a_vals[i]*bohr_to_angstrom:>10.4f} {c_vals[i]*bohr_to_angstrom:>10.4f} {energies[i]:>15.6f} {dE:>15.2f}")

print("-"*70)
print(f"Minimum energy: {E_min:.6f} eV/atom at")
print(f"  a = {a_vals[min_idx]*bohr_to_angstrom:.4f} Å ({a_vals[min_idx]:.4f} Bohr)")
print(f"  c = {c_vals[min_idx]*bohr_to_angstrom:.4f} Å ({c_vals[min_idx]:.4f} Bohr)")
