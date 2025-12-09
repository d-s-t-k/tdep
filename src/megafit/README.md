# MEGAFIT

## Short description

Fit force constants and thermodynamic properties across multi-dimensional parameter grids. MEGAFIT enables quasi-harmonic approximation (QHA) calculations with flexible parameter spaces including volume (V), temperature (T), and anisotropic lattice parameters (a, c). Force constants are interpolated using polynomial fitting, allowing computation of free energies at arbitrary points within the parameter space.

## Command line options

### Force constant cutoffs

| Option | Default | Description |
|--------|---------|-------------|
| `--secondorder_cutoff`, `-rc2` | `5.0` | Cutoff for the second order force constants (Å) |
| `--thirdorder_cutoff`, `-rc3` | `-1` | Cutoff for the third order force constants (Å). Negative means disabled. |
| `--fourthorder_cutoff`, `-rc4` | `-1` | Cutoff for the fourth order force constants (Å). Negative means disabled. |
| `--magnetic_pair_cutoff`, `-mc2` | `-1.0` | Cutoff for the pair magnetic interactions (Å) |

### Polar materials

| Option | Default | Description |
|--------|---------|-------------|
| `--polar` | `.false.` | Add dipole-dipole corrections for polar materials. Requires `infile.lotosplitting`. |
| `--polarcorrectiontype`, `-pc` | `3` | Type of polar correction to use. Choices: `1`, `2`, `3` |

### Fitting options

| Option | Default | Description |
|--------|---------|-------------|
| `--order`, `-o` | `2` | Order of the polynomials for the grid fitting procedure |
| `--pairfittype`, `-pf` | `1` | Method for fitting second order. `1` = global polynomial, `2` = locally adaptive polynomials |
| `--distancescale`, `-ds` | `0.1` | Scale factor for distances in the fitting procedure |
| `--temperaturescale`, `-ts` | `-1` | Scale factor for temperature. Negative uses automatic scaling. |

### Brillouin zone integration

| Option | Default | Description |
|--------|---------|-------------|
| `--harmonic_qpoint_grid`, `-qgh` | `26 26 26` | q-mesh density for harmonic (phonon) free energy calculations |
| `--anharmonic_qpoint_grid`, `-qga` | `10 10 10` | q-mesh density for anharmonic free energy calculations |

### Evaluation options

| Option | Default | Description |
|--------|---------|-------------|
| `--evalenergy` | `.false.` | Evaluate the free energy at points specified in `infile.evalpoints` |
| `--quasiharmonic` | `.false.` | In addition to the full anharmonic free energy, evaluate the quasiharmonic free energy for reference |

### General options

| Option | Default | Description |
|--------|---------|-------------|
| `--help`, `-h` | | Print help message |
| `--version`, `-v` | | Print version |

## Examples

```bash
# Basic QHA with 5 Å second-order cutoff
megafit -rc2 5.1

# Include third-order force constants
megafit -rc2 4.5 -rc3 3.21

# Polar material with dipole corrections
megafit -rc2 5.0 --polar

# Higher polynomial order for smoother interpolation
megafit -rc2 5.0 -o 3

# Finer q-mesh for converged phonon free energy
megafit -rc2 5.0 -qgh 30 30 30
```

## What does this code produce?

MEGAFIT interpolates force constants across a parameter grid (volume, temperature, lattice parameters) and computes thermodynamic quantities including:

- **Phonon free energy** $F_\text{ph}(x, T)$
- **Anharmonic corrections** from third and fourth order force constants
- **Total Helmholtz free energy** $F(x, T) = U(x) + \Delta U_0(x) + F_\text{ph}(x, T)$
- **Equilibrium parameters** (volume, lattice constants) as functions of temperature

The program supports several evaluation modes:

| Mode | Grid type | Description |
|------|-----------|-------------|
| 1 | Generic | User-defined grid evaluation |
| 3 | V-T | Volume-temperature grid with isobars |
| 4 | a-c | Lattice parameter grid (hexagonal/tetragonal) |
| 5 | a-c-T | Lattice parameters with temperature-dependent force constants |

## Theory

### Free energy components

The Helmholtz free energy is computed as:

$$F(x, T) = U(x) + \Delta U_0(x) + F_\text{ph}(x, T) + F_\text{ah}(x, T)$$

where:
- $U(x)$ = Static internal energy from equation of state
- $\Delta U_0(x)$ = Correction to static energy from force constant fitting
- $F_\text{ph}(x, T)$ = Phonon (harmonic) free energy
- $F_\text{ah}(x, T)$ = Anharmonic free energy corrections (if third/fourth order FCs provided)

### Phonon free energy

The phonon free energy per atom is:

$$F_\text{ph} = \frac{1}{N_a} \sum_{\mathbf{q},\nu} w_{\mathbf{q}} \left[ \frac{\hbar\omega_{\mathbf{q}\nu}}{2} + k_B T \ln\left(1 - e^{-\hbar\omega_{\mathbf{q}\nu}/k_B T}\right) \right]$$

where:
- $\omega_{\mathbf{q}\nu}$ = Phonon frequency at wavevector $\mathbf{q}$, branch $\nu$
- $w_{\mathbf{q}}$ = Integration weight
- $N_a$ = Number of atoms

### Polynomial interpolation

Force constants are interpolated across the parameter grid using polynomials:

$$\Phi_{ij}^{\alpha\beta}(x) = \sum_{n} C_{ij,n}^{\alpha\beta} P_n(x)$$

where $P_n(x)$ are polynomial basis functions and $x$ represents the grid coordinates (V, T, a, c, etc.).

For a-c grids, the free energy surface is fitted to a 4th order polynomial to find equilibrium:

$$F(a, c) = \sum_{i+j \leq 4} C_{ij} \tilde{a}^i \tilde{c}^j$$

## Input files

### Required files

| File | Description |
|------|-------------|
| `infile.ucposcar` | Unit cell structure in VASP POSCAR format |
| `infile.ssposcar` | Supercell structure in VASP POSCAR format |
| `infile.forceconstant` | Reference second-order force constants |
| `infile.simulations` | Parameter grid definition and paths to simulation data |
| `infile.meta` | Polynomial fit metadata (from `fitmultipole`) |

### Optional files

| File | Description |
|------|-------------|
| `infile.evalpoints` | Points at which to evaluate free energy |
| `infile.lotosplitting` | Born charges and dielectric tensor (polar materials) |
| `infile.forceconstant_thirdorder` | Third-order force constants |
| `infile.forceconstant_fourthorder` | Fourth-order force constants |

### `infile.simulations` format

This file defines the parameter grid and references the simulation data:

```
ndim                        # Number of dimensions (1, 2, or 3)
nsim                        # Number of simulation points
varnames                    # Variable names: V, T, a, c (space-separated)
order_per_dim               # Polynomial order per dimension
eosname                     # Equation of state: Birch, Vinet, 2D-Bi, or null
eos_params                  # EOS parameters (if using an EOS)
coord1 [static_E] /path/    # Grid coordinates, optional static energy (eV/atom), and path
coord2 [static_E] /path/
...
```

**Example for volume-only QHA (1D):**
```
1
5
V
3
Birch
-5.5 10.0 150.0 4.0
9.5  /path/to/sim1/outfile.grid_simulation.hdf5
9.8  /path/to/sim2/outfile.grid_simulation.hdf5
10.0 /path/to/sim3/outfile.grid_simulation.hdf5
10.2 /path/to/sim4/outfile.grid_simulation.hdf5
10.5 /path/to/sim5/outfile.grid_simulation.hdf5
```

**Example for a-c lattice parameter grid (2D) with static energies:**

For a-c grids with `eosname=null`, the static 0K DFT energy must be provided for each grid point.
The format includes an additional column for static energy in eV/atom:

```
2
9
a c
4 4
null
2.90 4.60 -127.456 /path/to/sim_a1_c1/outfile.sim.hdf5
2.90 4.65 -127.389 /path/to/sim_a1_c2/outfile.sim.hdf5
2.90 4.70 -127.312 /path/to/sim_a1_c3/outfile.sim.hdf5
2.95 4.60 -127.501 /path/to/sim_a2_c1/outfile.sim.hdf5
...
```

The static energy should be the total DFT energy of the relaxed primitive cell at each (a,c) point,
converted to eV per atom. For Quantum ESPRESSO, this is the "total energy" from the `ev/` calculation.
For VASP, this is the "energy without entropy" (E0) from OUTCAR.

### `infile.evalpoints` format

Defines where to evaluate the interpolated quantities:

**For volume-only QHA (evalmode=1):**
```
1                    # evalmode
nv lin               # Number of volume points, spacing type
Vmin Vmax            # Volume range (Å³/atom)
0                    # Pressure step (GPa), 0 = no pressure output
```

**For a-c lattice parameter QHA (evalmode=4):**
```
4                    # evalmode
na lin               # Number of a values, spacing
amin amax            # a range (Å)
nc lin               # Number of c values, spacing
cmin cmax            # c range (Å)
```

**For a-c-T grid with T-dependent FCs (evalmode=5):**
```
5                    # evalmode
na lin               # Number of a values, spacing
amin amax            # a range (Å)
nc lin               # Number of c values, spacing
cmin cmax            # c range (Å)
nt lin               # Number of T values, spacing
Tmin Tmax            # Temperature range (K)
```

## Output files

### `outfile.interpolated_free_energy.hdf5`

HDF5 file containing all computed thermodynamic quantities. The structure depends on the evaluation mode:

**For volume QHA (`grid_QHA` group):**
- `volumes` (Å³/atom)
- `temperatures` (K)
- `static_internal_energy` U(V) (eV/atom)
- `delta_U0` ΔU₀(V) (eV/atom)
- `phonon_free_energy` F_ph(V,T) (eV/atom)
- `Helmholtz_free_energy` F_total (eV/atom)

**For a-c QHA (`grid_ac_QHA` group):**
- `a_values`, `c_values` (Å)
- `temperatures` (K)
- `static_internal_energy` U(a,c) (eV/atom)
- `phonon_free_energy` F_ph(a,c,T) (eV/atom)
- `anharmonic_free_energy_3rd` (eV/atom)
- `anharmonic_free_energy_4th` (eV/atom)
- `Helmholtz_free_energy` F(a,c,T) (eV/atom)
- `a_equilibrium`, `c_equilibrium` (Å) - equilibrium vs T
- `F_equilibrium` (eV/atom) - minimum F at each T

### `outfile.ac_equilibrium.dat` / `outfile.act_equilibrium.dat`

Plain text file with equilibrium lattice parameters versus temperature:

```
# Temperature (K)    a_eq (Å)    c_eq (Å)    F_min (eV/atom)
1.000000E+00         2.950000E+00    4.650000E+00    -5.123456E+00
1.000000E+02         2.951234E+00    4.651234E+00    -5.122345E+00
...
```

## Workflow

### Step 1: Generate simulation data

Run DFT/MD simulations at multiple points on the parameter grid:

```bash
# For each grid point
cd sim_point_1/
# ... run simulation ...
extract_forceconstants -rc2 5.0
pack_simulation --output_format 1
cd ..
```

### Step 2: Fit force constants across the grid

```bash
fitmultipole
```

This produces `infile.meta` with polynomial coefficients.

### Step 3: Run megafit

```bash
megafit -rc2 5.0 -qgh 30 30 30
```

### Step 4: Analyze output

Use the HDF5 output for further analysis:

```python
import h5py
import numpy as np

with h5py.File('outfile.interpolated_free_energy.hdf5', 'r') as f:
    T = f['grid_ac_QHA/temperatures'][:]
    a_eq = f['grid_ac_QHA/a_equilibrium'][:]
    c_eq = f['grid_ac_QHA/c_equilibrium'][:]
    
    # Thermal expansion
    alpha_a = np.gradient(a_eq, T) / a_eq
    alpha_c = np.gradient(c_eq, T) / c_eq
```

## Tips and best practices

1. **Grid density**: Use at least 5 points per dimension for reliable polynomial interpolation

2. **q-mesh convergence**: Test convergence of phonon free energy with q-mesh density. Start with `26 26 26` and increase if needed.

3. **Polynomial order**: Higher orders capture more detail but may introduce oscillations. Order 2-3 is typically sufficient.

4. **Unstable modes**: Points with imaginary frequencies are flagged with F = 1.234×10⁸ and excluded from fitting

5. **Memory**: Large grids with fine q-meshes require significant memory. Use MPI for parallel execution:
   ```bash
   mpirun -np 8 megafit -rc2 5.0
   ```

6. **Debugging**: Use `--verbose` for detailed output during execution

## Troubleshooting

| Problem | Possible cause | Solution |
|---------|----------------|----------|
| "NOT DONE" error | evalmode mismatch | Check evalmode in `infile.evalpoints` matches grid type |
| Very large F values (10⁸) | Imaginary frequencies | Check structural stability at that grid point |
| Missing module files | Build order issue | Rebuild TDEP library: `./build_things.sh` |
| Memory error | Grid too large | Reduce q-mesh or use more MPI ranks |
| Poor interpolation | Insufficient grid | Add more simulation points |

## Related programs

- [`fitmultipole`](../fitmultipole/) - Fit force constants across parameter grid
- [`extract_forceconstants`](../extract_forceconstants/) - Extract force constants from MD
- [`phonon_dispersion_relations`](../phonon_dispersion_relations/) - Compute phonon dispersions  
- [`anharmonic_free_energy`](../anharmonic_free_energy/) - Compute anharmonic corrections
- [`pack_simulation`](../pack_simulation/) - Package simulation data to HDF5

## References

1. O. Hellman, P. Steneteg, I. A. Abrikosov, and S. I. Simak, "Temperature dependent effective potential method for accurate free energy calculations of solids", Phys. Rev. B **87**, 104111 (2013)

2. O. Hellman and I. A. Abrikosov, "Temperature-dependent effective third-order interatomic force constants from first principles", Phys. Rev. B **88**, 144301 (2013)

3. O. Hellman, I. A. Abrikosov, and S. I. Simak, "Lattice dynamics of anharmonic solids from first principles", Phys. Rev. B **84**, 180301(R) (2011)
