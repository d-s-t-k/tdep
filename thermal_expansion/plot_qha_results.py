#!/usr/bin/env python3
"""
Plot QHA results from megafit:
1. Interpolated vs training point free energies at selected temperatures
2. Lattice parameters (a, c) vs temperature
3. Free energy surface at different temperatures
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import RegularGridInterpolator

# Read the HDF5 output
h5_file = 'outfile.interpolated_free_energy.hdf5'

with h5py.File(h5_file, 'r') as f:
    grp = f['grid_ac_QHA']
    
    # Grid data
    a_values = grp['a_values'][:]  # in Angstrom
    c_values = grp['c_values'][:]  # in Angstrom
    temperatures = grp['temperatures'][:]  # in K
    
    # Free energy components on full grid 
    # HDF5 stores Fortran array F(na,nc,nt) as (nt,nc,na) in C order
    # Swap axes 1 and 2 to get (nt, na, nc) for Python indexing F[ti, ia, ic]
    F_total = np.swapaxes(grp['Helmholtz_free_energy'][:], 1, 2)
    F_phonon = np.swapaxes(grp['phonon_free_energy'][:], 1, 2)
    U_static = np.swapaxes(grp['static_internal_energy'][:], 1, 2)
    delta_U0 = np.swapaxes(grp['delta_U0'][:], 1, 2)
    
    # Equilibrium values (from polynomial fit - may not match grid minimum!)
    a_eq = grp['a_equilibrium'][:]
    c_eq = grp['c_equilibrium'][:]
    F_eq = grp['F_equilibrium'][:]

print(f"Grid: {len(a_values)} a x {len(c_values)} c = {len(a_values)*len(c_values)} points")
print(f"F_total shape after swapaxes: {F_total.shape} (nT, na, nc)")
print(f"a range: {a_values[0]:.4f} - {a_values[-1]:.4f} Å")
print(f"c range: {c_values[0]:.4f} - {c_values[-1]:.4f} Å")
print(f"Temperatures: {temperatures}")

# Parse training points from infile.simulations
training_a = []
training_c = []
training_E0 = []

with open('infile.simulations', 'r') as f:
    lines = f.readlines()
    for line in lines[5:]:  # Skip header lines
        parts = line.strip().split()
        if len(parts) >= 3:
            training_a.append(float(parts[0]))
            training_c.append(float(parts[1]))
            training_E0.append(float(parts[2]))

training_a = np.array(training_a)
training_c = np.array(training_c)
training_E0 = np.array(training_E0)

print(f"\nTraining points: {len(training_a)}")
print(f"Training a range: {training_a.min():.4f} - {training_a.max():.4f} Å")
print(f"Training c range: {training_c.min():.4f} - {training_c.max():.4f} Å")

# ============================================================
# Figure 1: Free energy comparison at different temperatures
# ============================================================
fig1, axes = plt.subplots(2, 2, figsize=(12, 10))

# Select temperatures to plot - pick 4 evenly spaced through the range
nT = len(temperatures)
temp_indices = [0, nT//3, 2*nT//3, nT-1]  # Low, mid-low, mid-high, high temperatures

for idx, (ax, ti) in enumerate(zip(axes.flat, temp_indices)):
    T = temperatures[ti]
    F_grid = F_total[ti, :, :]  # (na, nc)
    
    # Create interpolator for this temperature - allow bounds_error=False for edge points
    interp = RegularGridInterpolator((a_values, c_values), F_grid, 
                                      method='linear', bounds_error=False, fill_value=None)
    
    # Interpolate at training points
    F_interp_at_training = np.array([interp([a, c])[0] for a, c in zip(training_a, training_c)])
    
    # For comparison, we need to compute free energy at training points
    # The static energy at training points should match training_E0
    # Let's also get the static energy from interpolation
    U_static_grid = U_static[ti, :, :]
    interp_static = RegularGridInterpolator((a_values, c_values), U_static_grid, 
                                             method='linear', bounds_error=False, fill_value=None)
    U_static_at_training = np.array([interp_static([a, c])[0] for a, c in zip(training_a, training_c)])
    
    # Shift energies to have minimum at 0 for comparison
    F_min = F_interp_at_training.min()
    F_shifted = (F_interp_at_training - F_min) * 1000  # meV/atom
    
    # Create a 2D scatter plot with color by free energy
    sc = ax.scatter(training_a, training_c, c=F_shifted, cmap='viridis', 
                    s=100, edgecolor='black', linewidth=1)
    
    # Mark the equilibrium point (from polynomial fit)
    ax.scatter([a_eq[ti]], [c_eq[ti]], c='red', s=200, marker='*', 
               edgecolor='black', linewidth=1, zorder=5, label='Equilibrium')
    
    ax.set_xlabel('a (Å)', fontsize=11)
    ax.set_ylabel('c (Å)', fontsize=11)
    ax.set_title(f'T = {T:.0f} K', fontsize=12)
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('F - F_min (meV/atom)', fontsize=10)
    
    if idx == 0:
        ax.legend(loc='upper right')

plt.suptitle('Free Energy at Training Points\n(color = interpolated Helmholtz free energy)', fontsize=14)
plt.tight_layout()
plt.savefig('free_energy_at_training_points.png', dpi=150, bbox_inches='tight')
print("\nSaved: free_energy_at_training_points.png")

# ============================================================
# Figure 2: Free energy surface contours at different temperatures
# ============================================================
fig2, axes = plt.subplots(2, 2, figsize=(12, 10))

for idx, (ax, ti) in enumerate(zip(axes.flat, temp_indices)):
    T = temperatures[ti]
    F_grid = F_total[ti, :, :]
    F_min = F_grid.min()
    F_shifted = (F_grid - F_min) * 1000  # meV/atom
    
    # Create meshgrid for contour plot
    A, C = np.meshgrid(a_values, c_values, indexing='ij')
    
    # Contour plot
    levels = np.linspace(0, F_shifted.max() * 0.5, 20)  # Focus on low energy region
    cont = ax.contourf(A, C, F_shifted, levels=levels, cmap='viridis', extend='max')
    ax.contour(A, C, F_shifted, levels=levels, colors='white', linewidths=0.5, alpha=0.5)
    
    # Mark training points
    ax.scatter(training_a, training_c, c='white', s=40, edgecolor='black', 
               linewidth=0.5, marker='o', label='Training', zorder=4)
    
    # Mark equilibrium (from polynomial fit)
    ax.scatter([a_eq[ti]], [c_eq[ti]], c='red', s=200, marker='*', 
               edgecolor='white', linewidth=1.5, zorder=5, label='Equilibrium')
    
    ax.set_xlabel('a (Å)', fontsize=11)
    ax.set_ylabel('c (Å)', fontsize=11)
    ax.set_title(f'T = {T:.0f} K', fontsize=12)
    
    cbar = plt.colorbar(cont, ax=ax)
    cbar.set_label('F - F_min (meV/atom)', fontsize=10)
    
    if idx == 0:
        ax.legend(loc='upper right', fontsize=9)

plt.suptitle('Helmholtz Free Energy Surface', fontsize=14)
plt.tight_layout()
plt.savefig('free_energy_surface.png', dpi=150, bbox_inches='tight')
print("Saved: free_energy_surface.png")

# ============================================================
# Figure 3: Lattice parameters vs temperature
# ============================================================
fig3, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Plot a(T)
ax1 = axes[0]
ax1.plot(temperatures, a_eq, 'o-', color='blue', linewidth=2, markersize=6)
ax1.set_xlabel('Temperature (K)', fontsize=12)
ax1.set_ylabel('a (Å)', fontsize=12)
ax1.set_title('Lattice Parameter a vs Temperature', fontsize=12)
ax1.grid(True, alpha=0.3)

# Calculate thermal expansion coefficient for a
if len(temperatures) > 1:
    dT = temperatures[-1] - temperatures[0]
    da = a_eq[-1] - a_eq[0]
    a_avg = np.mean(a_eq)
    alpha_a = (da / a_avg) / dT * 1e6  # in 10^-6 K^-1
    ax1.text(0.05, 0.95, f'α_a ≈ {alpha_a:.1f} × 10⁻⁶ K⁻¹', 
             transform=ax1.transAxes, fontsize=10, verticalalignment='top')

# Plot c(T)
ax2 = axes[1]
ax2.plot(temperatures, c_eq, 's-', color='green', linewidth=2, markersize=6)
ax2.set_xlabel('Temperature (K)', fontsize=12)
ax2.set_ylabel('c (Å)', fontsize=12)
ax2.set_title('Lattice Parameter c vs Temperature', fontsize=12)
ax2.grid(True, alpha=0.3)

# Calculate thermal expansion coefficient for c
if len(temperatures) > 1:
    dc = c_eq[-1] - c_eq[0]
    c_avg = np.mean(c_eq)
    alpha_c = (dc / c_avg) / dT * 1e6  # in 10^-6 K^-1
    ax2.text(0.05, 0.95, f'α_c ≈ {alpha_c:.1f} × 10⁻⁶ K⁻¹', 
             transform=ax2.transAxes, fontsize=10, verticalalignment='top')

# Plot c/a ratio
ax3 = axes[2]
c_over_a = c_eq / a_eq
ax3.plot(temperatures, c_over_a, 'd-', color='purple', linewidth=2, markersize=8)
ax3.set_xlabel('Temperature (K)', fontsize=12)
ax3.set_ylabel('c/a ratio', fontsize=12)
ax3.set_title('c/a Ratio vs Temperature', fontsize=12)
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('lattice_parameters_vs_T.png', dpi=150, bbox_inches='tight')
print("Saved: lattice_parameters_vs_T.png")

# ============================================================
# Figure 4: Equilibrium free energy vs temperature
# ============================================================
fig4, ax = plt.subplots(figsize=(8, 6))

ax.plot(temperatures, F_eq * 1000, 'o-', color='crimson', linewidth=2, markersize=8)
ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('F_eq (meV/atom)', fontsize=12)
ax.set_title('Equilibrium Helmholtz Free Energy vs Temperature', fontsize=14)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('equilibrium_free_energy_vs_T.png', dpi=150, bbox_inches='tight')
print("Saved: equilibrium_free_energy_vs_T.png")

# ============================================================
# Figure 5: Interpolation quality check - 1D cuts through training points
# ============================================================
fig5, axes = plt.subplots(2, 3, figsize=(14, 8))

# Get unique a and c values from training set
unique_a = np.unique(training_a)
unique_c = np.unique(training_c)

# Plot 1D cuts at constant c (varying a) for different temperatures
# Pick 3 evenly spaced temperatures
cut_temp_indices = [0, nT//2, nT-1]
for ti, ax in zip(cut_temp_indices, axes[0, :]):
    T = temperatures[ti]
    F_grid = F_total[ti, :, :]
    
    # Use middle c value from training
    c_mid = unique_c[len(unique_c)//2]
    c_idx = np.argmin(np.abs(c_values - c_mid))
    
    # Get the 1D cut
    F_cut = F_grid[:, c_idx]
    F_min = F_cut.min()
    
    ax.plot(a_values, (F_cut - F_min) * 1000, 'b-', linewidth=2, label='Interpolated')
    
    # Training points at this c
    mask = np.abs(training_c - c_mid) < 0.01
    a_train = training_a[mask]
    
    # Get interpolated values at training a for comparison
    interp = RegularGridInterpolator((a_values, c_values), F_grid, 
                                      method='linear', bounds_error=False, fill_value=None)
    F_at_train = np.array([interp([a, c_mid])[0] for a in a_train])
    
    ax.scatter(a_train, (F_at_train - F_min) * 1000, c='red', s=80, 
               edgecolor='black', zorder=5, label='Training points')
    
    ax.set_xlabel('a (Å)', fontsize=11)
    ax.set_ylabel('F - F_min (meV/atom)', fontsize=11)
    ax.set_title(f'T = {T:.0f} K, c = {c_mid:.4f} Å', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

# Plot 1D cuts at constant a (varying c) for different temperatures
for ti, ax in zip(cut_temp_indices, axes[1, :]):
    T = temperatures[ti]
    F_grid = F_total[ti, :, :]
    
    # Use middle a value from training
    a_mid = unique_a[len(unique_a)//2]
    a_idx = np.argmin(np.abs(a_values - a_mid))
    
    # Get the 1D cut
    F_cut = F_grid[a_idx, :]
    F_min = F_cut.min()
    
    ax.plot(c_values, (F_cut - F_min) * 1000, 'g-', linewidth=2, label='Interpolated')
    
    # Training points at this a
    mask = np.abs(training_a - a_mid) < 0.01
    c_train = training_c[mask]
    
    # Get interpolated values at training c for comparison
    interp = RegularGridInterpolator((a_values, c_values), F_grid, 
                                      method='linear', bounds_error=False, fill_value=None)
    F_at_train = np.array([interp([a_mid, c])[0] for c in c_train])
    
    ax.scatter(c_train, (F_at_train - F_min) * 1000, c='red', s=80, 
               edgecolor='black', zorder=5, label='Training points')
    
    ax.set_xlabel('c (Å)', fontsize=11)
    ax.set_ylabel('F - F_min (meV/atom)', fontsize=11)
    ax.set_title(f'T = {T:.0f} K, a = {a_mid:.4f} Å', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

plt.suptitle('1D Free Energy Cuts Through Training Points', fontsize=14)
plt.tight_layout()
plt.savefig('interpolation_1d_cuts.png', dpi=150, bbox_inches='tight')
print("Saved: interpolation_1d_cuts.png")

# ============================================================
# Print summary
# ============================================================
print("\n" + "="*60)
print("SUMMARY OF QHA RESULTS")
print("="*60)
print(f"\nTemperature range: {temperatures[0]:.0f} - {temperatures[-1]:.0f} K")
print(f"\nEquilibrium lattice parameters:")
print(f"  T = {temperatures[0]:.0f} K:  a = {a_eq[0]:.5f} Å,  c = {c_eq[0]:.5f} Å,  c/a = {c_eq[0]/a_eq[0]:.4f}")
print(f"  T = {temperatures[-1]:.0f} K:  a = {a_eq[-1]:.5f} Å,  c = {c_eq[-1]:.5f} Å,  c/a = {c_eq[-1]/a_eq[-1]:.4f}")
print(f"\nChange from {temperatures[0]:.0f} K to {temperatures[-1]:.0f} K:")
print(f"  Δa = {(a_eq[-1] - a_eq[0])*1000:.3f} mÅ  ({(a_eq[-1]/a_eq[0] - 1)*100:.3f}%)")
print(f"  Δc = {(c_eq[-1] - c_eq[0])*1000:.3f} mÅ  ({(c_eq[-1]/c_eq[0] - 1)*100:.3f}%)")
print(f"\nLinear thermal expansion coefficients (averaged):")
print(f"  α_a ≈ {alpha_a:.1f} × 10⁻⁶ K⁻¹")
print(f"  α_c ≈ {alpha_c:.1f} × 10⁻⁶ K⁻¹")
print(f"\nEquilibrium free energies:")
print(f"  F({temperatures[0]:.0f} K) = {F_eq[0]*1000:.3f} meV/atom")
print(f"  F({temperatures[-1]:.0f} K) = {F_eq[-1]*1000:.3f} meV/atom")

plt.show()
