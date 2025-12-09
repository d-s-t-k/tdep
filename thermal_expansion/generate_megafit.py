#!/usr/bin/env python3
"""
Generate infile.simulations and infile.evalpoints for megafit QHA.
Scans a*_c* folders, extracts lattice parameters and energies.
"""

import os
import re
import subprocess
import numpy as np
from pathlib import Path

# Constants
RY_TO_EV = 13.605693122994  # Rydberg to eV

def get_lattice_params_from_ucposcar(ucposcar_path):
    """
    Extract hexagonal a and c from rhombohedral ucposcar.
    For rhombohedral cell with vectors v1, v2, v3:
    - a_hex = |v1 - v2| (in-plane)
    - c_hex = |v1 + v2 + v3| (along [111])
    """
    with open(ucposcar_path, 'r') as f:
        lines = f.readlines()
    
    scale = float(lines[1].strip())
    v1 = np.array([float(x) for x in lines[2].split()]) * scale
    v2 = np.array([float(x) for x in lines[3].split()]) * scale
    v3 = np.array([float(x) for x in lines[4].split()]) * scale
    
    # Hexagonal a = distance between equivalent atoms in basal plane
    a_hex = np.linalg.norm(v1 - v2)
    # Hexagonal c = 3 * height along [111] direction  
    c_hex = np.linalg.norm(v1 + v2 + v3)
    
    return a_hex, c_hex

def get_energy_from_qe(ev_folder):
    """Extract total energy from Quantum Espresso output."""
    pwo_file = os.path.join(ev_folder, 'espresso.pwo')
    if not os.path.exists(pwo_file):
        return None
    
    energy_ry = None
    with open(pwo_file, 'r') as f:
        for line in f:
            if '!    total energy' in line:
                # Extract energy in Ry
                match = re.search(r'=\s+([-\d.]+)\s+Ry', line)
                if match:
                    energy_ry = float(match.group(1))
    
    if energy_ry is not None:
        return energy_ry * RY_TO_EV
    return None

def get_natoms_from_ucposcar(ucposcar_path):
    """Get number of atoms from ucposcar."""
    with open(ucposcar_path, 'r') as f:
        lines = f.readlines()
    # Line 7 has atom counts
    counts = [int(x) for x in lines[6].split()]
    return sum(counts)

def run_pack_simulation(folder):
    """Run pack_simulation in the given folder."""
    hdf5_file = os.path.join(folder, 'outfile.sim.hdf5')
    if os.path.exists(hdf5_file):
        print(f"  {folder}: outfile.sim.hdf5 already exists, skipping pack_simulation")
        return True
    
    print(f"  {folder}: running pack_simulation...")
    try:
        result = subprocess.run(
            ['pack_simulation', '-of', 'hdf5'],
            cwd=folder,
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode != 0:
            print(f"    WARNING: pack_simulation failed: {result.stderr}")
            return False
        return True
    except Exception as e:
        print(f"    ERROR: {e}")
        return False

def main():
    base_dir = Path(__file__).parent.resolve()
    
    # Find all a*_c* folders
    folders = sorted([d for d in os.listdir(base_dir) 
                      if os.path.isdir(os.path.join(base_dir, d)) 
                      and re.match(r'a[\d.]+_c[\d.]+', d)])
    
    if not folders:
        print("ERROR: No a*_c* folders found!")
        return
    
    print(f"Found {len(folders)} simulation folders")
    
    # Collect simulation data
    simulations = []
    a_values = set()
    c_values = set()
    
    for folder in folders:
        folder_path = os.path.join(base_dir, folder)
        ucposcar = os.path.join(folder_path, 'infile.ucposcar')
        ev_folder = os.path.join(folder_path, 'ev')
        
        if not os.path.exists(ucposcar):
            print(f"  WARNING: {folder} missing infile.ucposcar, skipping")
            continue
        
        # Run pack_simulation
        if not run_pack_simulation(folder_path):
            continue
        
        # Get lattice parameters
        a, c = get_lattice_params_from_ucposcar(ucposcar)
        a_values.add(a)
        c_values.add(c)
        
        # Get energy
        energy_ev = get_energy_from_qe(ev_folder)
        if energy_ev is None:
            print(f"  WARNING: {folder} no energy found in ev/espresso.pwo")
            energy_ev = 0.0  # Will need manual correction
        
        # Get natoms for eV/atom
        natoms = get_natoms_from_ucposcar(ucposcar)
        energy_per_atom = energy_ev / natoms
        
        hdf5_path = os.path.join(folder_path, 'outfile.sim.hdf5')
        
        simulations.append({
            'a': a,
            'c': c,
            'energy': energy_per_atom,
            'hdf5': hdf5_path
        })
        
        print(f"  {folder}: a={a:.6f}, c={c:.6f}, E={energy_per_atom:.8f} eV/atom")
    
    if not simulations:
        print("ERROR: No valid simulations found!")
        return
    
    # Sort by a, then c
    simulations.sort(key=lambda x: (x['a'], x['c']))
    
    # Determine grid dimensions
    a_sorted = sorted(a_values)
    c_sorted = sorted(c_values)
    na = len(a_sorted)
    nc = len(c_sorted)
    
    print(f"\nGrid: {na} x {nc} (a x c)")
    print(f"a range: {min(a_sorted):.6f} to {max(a_sorted):.6f} Bohr")
    print(f"c range: {min(c_sorted):.6f} to {max(c_sorted):.6f} Bohr")
    
    # Write infile.simulations
    sim_file = os.path.join(base_dir, 'infile.simulations')
    with open(sim_file, 'w') as f:
        f.write("2\n")  # number of dimensions (a, c)
        f.write(f"{len(simulations)}\n")  # number of simulations
        f.write("a c\n")  # dimension names
        f.write("2 2\n")  # polynomial order per dimension
        f.write("null\n")  # reference structure (use first simulation)
        for sim in simulations:
            f.write(f"{sim['a']:.6f} {sim['c']:.6f} {sim['energy']:.8f} \"{sim['hdf5']}\"\n")
    
    print(f"\nWrote {sim_file}")
    
    # Write infile.evalpoints (for a-c grid mode)
    # Add some padding to the grid range
    a_min, a_max = min(a_sorted), max(a_sorted)
    c_min, c_max = min(c_sorted), max(c_sorted)
    a_pad = (a_max - a_min) * 0.05
    c_pad = (c_max - c_min) * 0.05
    
    eval_file = os.path.join(base_dir, 'infile.evalpoints')
    with open(eval_file, 'w') as f:
        f.write("4\n")  # evalmode 4 = a-c grid
        f.write("21 linear\n")  # 21 points for a
        f.write(f"{a_min - a_pad:.2f} {a_max + a_pad:.2f}\n")
        f.write("21 linear\n")  # 21 points for c
        f.write(f"{c_min - c_pad:.2f} {c_max + c_pad:.2f}\n")
    
    print(f"Wrote {eval_file}")
    
    # Print the megafit command
    print("\n" + "="*60)
    print("Run megafit with:")
    print("="*60)
    print(f"\ncd {base_dir}")
    print("megafit --evalenergy --quasiharmonic --temperature_range 1 1000 40 -qgh 8 8 8")
    print("\nOptions:")
    print("  --temperature_range MIN MAX NPTS  (e.g., 1 1000 40 for 1-1000K)")
    print("  -qgh NX NY NZ                     (harmonic q-point grid for phonons)")
    print()

if __name__ == '__main__':
    main()
