#!/usr/bin/env python3
"""
Compare interpolated free energies vs. training point free energies.
This validates the polynomial interpolation quality in megafit.
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
from pathlib import Path

# Paths
base_dir = Path('/Users/ds.kim/Software/tdep/thermal_expansion')
hdf5_file = base_dir / 'outfile.interpolated_free_energy.hdf5'
sim_file = base_dir / 'infile.simulations'

print("="*70)
print("INTERPOLATION QUALITY ANALYSIS")
print("="*70)

#------------------------------------------------------------------------------
# 1. Read training point data from infile.simulations
#------------------------------------------------------------------------------
print("\n1. Reading training point data...")

with open(sim_file, 'r') as f:
    lines = f.readlines()

ndim = int(lines[0].strip())
nsim = int(lines[1].strip())

# Parse simulation data
train_a = []
train_c = []
train_E_static = []  # Static DFT energy (eV/atom)
train_paths = []

for i in range(5, 5 + nsim):
    parts = lines[i].split('"')
    coords = parts[0].split()
    train_a.append(float(coords[0]))
    train_c.append(float(coords[1]))
    train_E_static.append(float(coords[2]))
    train_paths.append(parts[1])

train_a = np.array(train_a)
train_c = np.array(train_c)
train_E_static = np.array(train_E_static)

print(f"  Found {nsim} training points")
print(f"  a range: {train_a.min():.4f} to {train_a.max():.4f} Å")
print(f"  c range: {train_c.min():.4f} to {train_c.max():.4f} Å")

#------------------------------------------------------------------------------
# 2. Read interpolated data from HDF5
#------------------------------------------------------------------------------
print("\n2. Reading interpolated data from HDF5...")

with h5py.File(hdf5_file, 'r') as h5:
    grp = h5['grid_ac_QHA']
    
    # Grid axes
    a_grid = grp['a_values'][:]
    c_grid = grp['c_values'][:]
    temperatures = grp['temperatures'][:]
    
    # Energy components (T, a, c) in eV/atom
    F_total = grp['Helmholtz_free_energy'][:]    # (nT, na, nc)
    F_phonon = grp['phonon_free_energy'][:]
    U_static = grp['static_internal_energy'][:]
    U0 = grp['delta_U0'][:]
    ah3 = grp['anharmonic_free_energy_3rd'][:]
    ah4 = grp['anharmonic_free_energy_4th'][:]
    
    # Polynomial coefficients
    poly_coeffs = grp['polynomial_coefficients'][:]  # (nT, 25)

nT, na, nc = F_total.shape
print(f"  Grid: {na} a × {nc} c × {nT} T")
print(f"  a range: {a_grid.min():.4f} to {a_grid.max():.4f} Å")
print(f"  c range: {c_grid.min():.4f} to {c_grid.max():.4f} Å")
print(f"  T range: {temperatures.min():.1f} to {temperatures.max():.1f} K")

#------------------------------------------------------------------------------
# 3. Find interpolated values at training points
#------------------------------------------------------------------------------
print("\n3. Comparing interpolated vs training static energies...")

def find_nearest_indices(a_val, c_val, a_grid, c_grid):
    """Find nearest grid indices for a given (a, c) point."""
    ia = np.argmin(np.abs(a_grid - a_val))
    ic = np.argmin(np.abs(c_grid - c_val))
    return ia, ic

def bilinear_interp(a_val, c_val, a_grid, c_grid, data_2d):
    """Bilinear interpolation on 2D grid."""
    # Find bracketing indices
    ia = np.searchsorted(a_grid, a_val) - 1
    ic = np.searchsorted(c_grid, c_val) - 1
    
    # Clamp to valid range
    ia = np.clip(ia, 0, len(a_grid) - 2)
    ic = np.clip(ic, 0, len(c_grid) - 2)
    
    # Interpolation weights
    a0, a1 = a_grid[ia], a_grid[ia + 1]
    c0, c1 = c_grid[ic], c_grid[ic + 1]
    
    wa = (a_val - a0) / (a1 - a0) if a1 != a0 else 0.0
    wc = (c_val - c0) / (c1 - c0) if c1 != c0 else 0.0
    
    # Bilinear interpolation
    v00 = data_2d[ia, ic]
    v01 = data_2d[ia, ic + 1]
    v10 = data_2d[ia + 1, ic]
    v11 = data_2d[ia + 1, ic + 1]
    
    return (1-wa)*(1-wc)*v00 + (1-wa)*wc*v01 + wa*(1-wc)*v10 + wa*wc*v11

# Compare static energies at training points
interp_E_static = []
for i in range(nsim):
    # Use first temperature slice (static energy is T-independent)
    E_interp = bilinear_interp(train_a[i], train_c[i], a_grid, c_grid, U_static[0, :, :])
    interp_E_static.append(E_interp)

interp_E_static = np.array(interp_E_static)
error_static = interp_E_static - train_E_static

print(f"\n  Static Energy Comparison (interpolated - training):")
print(f"    Mean error:     {np.mean(error_static)*1000:.4f} meV/atom")
print(f"    Max error:      {np.max(np.abs(error_static))*1000:.4f} meV/atom")
print(f"    RMS error:      {np.sqrt(np.mean(error_static**2))*1000:.4f} meV/atom")

#------------------------------------------------------------------------------
# 4. Check phonon free energies at training points (if available)
#------------------------------------------------------------------------------
print("\n4. Checking phonon free energies at training point folders...")

# Try to read phonon free energies from training folders
train_F_phonon = {}  # {folder: {T: F_ph}}
sample_temps = [300, 500, 800]  # K

for i, path in enumerate(train_paths):
    folder = Path(path).parent
    free_energy_file = folder / 'outfile.free_energy'
    
    if free_energy_file.exists():
        with open(free_energy_file, 'r') as f:
            content = f.read().strip()
        # Parse: typically "T F_ph" format
        # Format varies, try to read
        try:
            lines_fe = content.split('\n')
            for line in lines_fe:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        T_train = float(parts[0])
                        F_train = float(parts[1])
                        if i not in train_F_phonon:
                            train_F_phonon[i] = {}
                        train_F_phonon[i][T_train] = F_train
        except:
            pass

if train_F_phonon:
    print(f"  Found phonon free energy data in {len(train_F_phonon)} training folders")
else:
    print("  No outfile.free_energy found in training folders")

#------------------------------------------------------------------------------
# 5. Create comparison plots
#------------------------------------------------------------------------------
print("\n5. Creating comparison plots...")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Plot 1: Static energy surface (interpolated)
ax = axes[0, 0]
A, C = np.meshgrid(a_grid, c_grid, indexing='ij')
contour = ax.contourf(A, C, U_static[0, :, :], levels=30, cmap='viridis')
ax.scatter(train_a, train_c, c='red', s=50, edgecolors='white', label='Training points')
ax.set_xlabel('a (Å)')
ax.set_ylabel('c (Å)')
ax.set_title('Interpolated Static Energy U(a,c)')
plt.colorbar(contour, ax=ax, label='U (eV/atom)')

# Plot 2: Training vs interpolated static energy
ax = axes[0, 1]
ax.scatter(train_E_static, interp_E_static, c='blue', s=50, alpha=0.7)
E_min, E_max = train_E_static.min(), train_E_static.max()
ax.plot([E_min, E_max], [E_min, E_max], 'r--', linewidth=2, label='Perfect fit')
ax.set_xlabel('Training Static Energy (eV/atom)')
ax.set_ylabel('Interpolated Static Energy (eV/atom)')
ax.set_title('Static Energy: Training vs Interpolated')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Error distribution
ax = axes[0, 2]
ax.hist(error_static * 1000, bins=15, color='steelblue', edgecolor='black')
ax.axvline(x=0, color='red', linestyle='--', linewidth=2)
ax.set_xlabel('Error (meV/atom)')
ax.set_ylabel('Count')
ax.set_title(f'Static Energy Error Distribution\nRMS = {np.sqrt(np.mean(error_static**2))*1000:.3f} meV/atom')

# Plot 4: Phonon free energy surface at T=300K
ax = axes[1, 0]
iT_300 = np.argmin(np.abs(temperatures - 300))
contour = ax.contourf(A, C, F_phonon[iT_300, :, :], levels=30, cmap='plasma')
ax.scatter(train_a, train_c, c='white', s=50, edgecolors='black', label='Training points')
ax.set_xlabel('a (Å)')
ax.set_ylabel('c (Å)')
ax.set_title(f'Phonon Free Energy F_ph(a,c) at T={temperatures[iT_300]:.0f}K')
plt.colorbar(contour, ax=ax, label='F_ph (eV/atom)')

# Plot 5: Total free energy surface at T=300K
ax = axes[1, 1]
contour = ax.contourf(A, C, F_total[iT_300, :, :], levels=30, cmap='RdYlBu_r')
ax.scatter(train_a, train_c, c='black', s=50, edgecolors='white', label='Training points')

# Mark minimum
F_slice = F_total[iT_300, :, :]
min_idx = np.unravel_index(np.argmin(F_slice), F_slice.shape)
ax.scatter(a_grid[min_idx[0]], c_grid[min_idx[1]], c='lime', s=200, marker='*', 
          edgecolors='black', linewidths=2, label='Minimum')
ax.set_xlabel('a (Å)')
ax.set_ylabel('c (Å)')
ax.set_title(f'Total Free Energy F(a,c) at T={temperatures[iT_300]:.0f}K')
plt.colorbar(contour, ax=ax, label='F (eV/atom)')
ax.legend()

# Plot 6: Free energy vs temperature at equilibrium
ax = axes[1, 2]
with h5py.File(hdf5_file, 'r') as h5:
    grp = h5['grid_ac_QHA']
    a_eq = grp['a_equilibrium'][:]
    c_eq = grp['c_equilibrium'][:]
    F_eq = grp['F_equilibrium'][:]

ax.plot(temperatures, F_eq, 'b-', linewidth=2)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Free Energy at Equilibrium (eV/atom)')
ax.set_title('F(T) at Equilibrium (a*, c*)')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(base_dir / 'interpolation_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"\nSaved: interpolation_comparison.png")

#------------------------------------------------------------------------------
# 6. Detailed comparison table
#------------------------------------------------------------------------------
print("\n" + "="*70)
print("DETAILED TRAINING POINT COMPARISON")
print("="*70)
print(f"{'#':>3} {'a (Å)':>10} {'c (Å)':>10} {'U_train':>14} {'U_interp':>14} {'Error (meV)':>12}")
print("-"*70)

for i in range(nsim):
    err = error_static[i] * 1000
    print(f"{i+1:>3} {train_a[i]:>10.4f} {train_c[i]:>10.4f} {train_E_static[i]:>14.6f} {interp_E_static[i]:>14.6f} {err:>12.3f}")

print("-"*70)
print(f"{'':>3} {'':>10} {'':>10} {'':>14} {'RMS Error:':>14} {np.sqrt(np.mean(error_static**2))*1000:>12.3f}")
