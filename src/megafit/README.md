# MEGAFIT - Thermodynamic Property Interpolation

## Overview

MEGAFIT is a program for computing thermodynamic properties (free energies) by interpolating force constants across parameter grids. It enables quasi-harmonic approximation (QHA) calculations with flexible parameter spaces including volume, temperature, and lattice parameters (a, c).

## Prerequisites

Before running `megafit`, you need:

1. **Computed force constants at grid points**: Run DFT/MD simulations at various parameter combinations and extract force constants using `extract_forceconstants`
2. **Fitted polynomial coefficients**: Run `fitmultipole` to fit the force constants across the parameter grid
3. **Input files** (described below)

## Input Files

### 1. `infile.simulations` (required)

Defines the parameter grid, equation of state, and references simulation data.

**Full format:**
```
ndim                        # Number of dimensions (1, 2, or 3)
nsim                        # Number of simulation points
varnames                    # Variable names (V, T, eta, a, c) space-separated
order_per_dim               # Polynomial order per dimension (e.g., "3 3" for 2D)
eosname                     # Equation of state: Birch, Vinet, 2D-Bi, or null
eos_params                  # EOS parameters (see below)
coord1 /path/to/sim1/       # Grid coordinates and path to simulation HDF5 file
coord2 /path/to/sim2/
...
```

**Equation of State (EOS) options:**
- `Birch`: Birch-Murnaghan EOS (4 parameters: E0, V0, B0, B0')
- `Vinet`: Vinet EOS (4 parameters: E0, V0, B0, B0')
- `2D-Bi`: 2D Birch-Murnaghan (9 parameters)
- `null`: No EOS (static energy set to zero)

The EOS provides the **static internal energy** U(V) or U(V,eta).

**Example for volume-only QHA (1D grid):**
```
1                           # Number of dimensions
5                           # Number of simulations
V                           # Dimension name (V = volume)
3                           # Polynomial order for V
Birch                       # Birch-Murnaghan EOS
-5.5 10.0 150.0 4.0         # E0 (eV), V0 (Å³), B0 (GPa), B0'
9.5  /path/to/sim1/outfile.grid_simulation.hdf5
9.8  /path/to/sim2/outfile.grid_simulation.hdf5
10.0 /path/to/sim3/outfile.grid_simulation.hdf5
10.2 /path/to/sim4/outfile.grid_simulation.hdf5
10.5 /path/to/sim5/outfile.grid_simulation.hdf5
```

**Example for a-c lattice parameter grid (2D grid):**
```
2                           # Number of dimensions  
9                           # Number of simulations
a c                         # Dimension names (lattice parameters in Angstrom)
3 3                         # Polynomial order for a, c
null                        # No EOS (use interpolated energy)
2.90 4.60 /path/to/sim_a1_c1/outfile.grid_simulation.hdf5
2.90 4.65 /path/to/sim_a1_c2/outfile.grid_simulation.hdf5
2.90 4.70 /path/to/sim_a1_c3/outfile.grid_simulation.hdf5
2.95 4.60 /path/to/sim_a2_c1/outfile.grid_simulation.hdf5
...
```

**Example for a-c-T grid with T-dependent FCs (3D grid):**
```
3                           # Number of dimensions  
27                          # Number of simulations (3×3×3)
a c T                       # Dimension names
3 3 2                       # Polynomial order for a, c, T
null                        # No EOS
2.90 4.60 100 /path/to/sim_a1_c1_T1/outfile.grid_simulation.hdf5
2.90 4.60 300 /path/to/sim_a1_c1_T2/outfile.grid_simulation.hdf5
2.90 4.60 600 /path/to/sim_a1_c1_T3/outfile.grid_simulation.hdf5
2.90 4.65 100 /path/to/sim_a1_c2_T1/outfile.grid_simulation.hdf5
...
```

### 2. `infile.evalpoints` (required)

Defines where to evaluate the interpolated free energy.

**For volume-only QHA (evalmode=1):**
```
1                    # evalmode (1 = generic grid definition)
nv lin               # Number of volume points, spacing type (lin/den)
Vmin Vmax            # Volume range in Å³/atom
0                    # Pressure step (GPa), 0 = no pressure output
```

**For a-c lattice parameter QHA (evalmode=4):**
FCs interpolated in (a,c), free energy evaluated at many temperatures.
```
4                    # evalmode (4 = a-c lattice parameter grid)
na lin               # Number of a values, spacing type
amin amax            # a range in Angstrom
nc lin               # Number of c values, spacing type
cmin cmax            # c range in Angstrom
```

**For a-c-T grid with T-dependent FCs (evalmode=5):**
FCs interpolated in (a,c,T), free energy evaluated at each (a,c,T) point.
```
5                    # evalmode (5 = a-c-T grid with T-dependent FCs)
na lin               # Number of a values, spacing type
amin amax            # a range in Angstrom
nc lin               # Number of c values, spacing type
cmin cmax            # c range in Angstrom
nt lin               # Number of temperature values, spacing type
Tmin Tmax            # Temperature range in K
```

**For V-T grid (evalmode=3):**
```
3                    # evalmode (3 = V-T grid with adaptive volume)
nt lin               # Number of temperature points, spacing type
Tmin Tmax            # Temperature range in K
nv lin               # Number of volume points, spacing type
Pstep                # Pressure step in GPa
```

### 3. `infile.ucposcar` (required)

Reference crystal structure in VASP POSCAR format.

### 4. `infile.forceconstant` (required)

Reference second-order force constants (from `extract_forceconstants`).

### 5. `infile.meta` (required)

Metadata file containing polynomial fit coefficients (from `fitmultipole`).

### 6. `infile.forces` (optional)

Reference forces for higher-order corrections.

### 7. `infile.lotosplitting` (optional)

LO-TO splitting data for polar materials.

## Command Line Options

```bash
megafit [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--qmesh_density` | Q-mesh density for harmonic calculations | 26 26 26 |
| `--readforcemap` | Read forcemap from file | false |
| `--dumpgrid` | Dump force constants at all grid points | false |
| `--help`, `-h` | Show help message | |

## Usage Examples

### Example 1: Volume-only QHA (isotropic materials)

1. **Run DFT calculations** at 5+ volumes (e.g., -4%, -2%, 0%, +2%, +4%)

2. **Extract force constants** for each volume:
   ```bash
   cd sim_V1/ && extract_forceconstants && cd ..
   # repeat for each volume
   ```

3. **Create `infile.simulations`**:
   ```
   1
   V
   5
   9.5  ./sim_V1/
   9.8  ./sim_V2/
   10.0 ./sim_V3/
   10.2 ./sim_V4/
   10.5 ./sim_V5/
   ```

4. **Run fitmultipole** to fit force constants:
   ```bash
   fitmultipole
   ```

5. **Create `infile.evalpoints`**:
   ```
   1
   50 lin
   9.0 11.0
   0
   ```

6. **Run megafit**:
   ```bash
   megafit --qmesh_density 30 30 30
   ```

### Example 2: a-c Lattice Parameter QHA (hexagonal/tetragonal)

1. **Run DFT calculations** on a grid of (a, c) values:
   ```
   (a1,c1) (a1,c2) (a1,c3)
   (a2,c1) (a2,c2) (a2,c3)
   (a3,c1) (a3,c2) (a3,c3)
   ```

2. **Extract force constants** for each (a,c) point

3. **Create `infile.simulations`**:
   ```
   2
   a c
   9
   2.90 4.60 ./sim_a1_c1/
   2.90 4.65 ./sim_a1_c2/
   2.90 4.70 ./sim_a1_c3/
   2.95 4.60 ./sim_a2_c1/
   2.95 4.65 ./sim_a2_c2/
   2.95 4.70 ./sim_a2_c3/
   3.00 4.60 ./sim_a3_c1/
   3.00 4.65 ./sim_a3_c2/
   3.00 4.70 ./sim_a3_c3/
   ```

4. **Run fitmultipole** to fit force constants

5. **Create `infile.evalpoints`**:
   ```
   4
   25 lin
   2.85 3.05
   25 lin
   4.55 4.75
   ```

6. **Run megafit**:
   ```bash
   megafit --qmesh_density 30 30 30
   ```

### Example 3: V-T Grid (temperature-dependent volumes)

1. **Run MD calculations** at multiple (V,T) combinations

2. **Create `infile.simulations`**:
   ```
   2
   V T
   20
   9.5 100 ./sim_V1_T1/
   9.5 300 ./sim_V1_T2/
   ...
   ```

3. **Create `infile.evalpoints`**:
   ```
   1
   100 lin
   10 2000
   50 lin
   9.0 11.0
   0.5
   ```

4. **Run megafit**:
   ```bash
   megafit --qmesh_density 30 30 30
   ```

## Output Files

### `outfile.interpolated_free_energy.hdf5`

HDF5 file containing:

**For volume QHA (`grid_QHA` group):**
- `volumes`: Volume axis (Å³/atom)
- `temperatures`: Temperature axis (K)
- `static_internal_energy`: U(V) from EOS (eV/atom)
- `delta_U0`: Correction to static energy (eV/atom)
- `phonon_free_energy`: F_ph(V,T) (eV/atom)
- `Helmholtz_free_energy`: F_total = U + ΔU₀ + F_ph (eV/atom)

**For a-c QHA (`grid_ac_QHA` group):**
- `a_values`: a lattice parameter axis (Å)
- `c_values`: c lattice parameter axis (Å)
- `temperatures`: Temperature axis (K)
- `static_internal_energy`: U(a,c) (eV/atom)
- `phonon_free_energy`: F_ph(a,c,T) (eV/atom)
- `anharmonic_free_energy_3rd`: F_ah3(a,c,T) (eV/atom)
- `anharmonic_free_energy_4th`: F_ah4(a,c,T) (eV/atom)
- `Helmholtz_free_energy`: F_total(a,c,T) (eV/atom)
- `a_equilibrium`: Equilibrium a(T) (Å)
- `c_equilibrium`: Equilibrium c(T) (Å)
- `F_equilibrium`: Minimum F at each T (eV/atom)
- `polynomial_coefficients`: 4th order fit coefficients

**For a-c-T with T-dependent FCs (`grid_act` group):**
- `a_values`: a lattice parameter axis (Å)
- `c_values`: c lattice parameter axis (Å)
- `temperatures`: Temperature axis (K)
- `static_internal_energy`: U(a,c,T) (eV/atom)
- `phonon_free_energy`: F_ph(a,c,T) (eV/atom)
- `anharmonic_free_energy_3rd`: F_ah3(a,c,T) (eV/atom)
- `anharmonic_free_energy_4th`: F_ah4(a,c,T) (eV/atom)
- `Helmholtz_free_energy`: F_total(a,c,T) (eV/atom)
- `a_equilibrium`: Equilibrium a(T) (Å)
- `c_equilibrium`: Equilibrium c(T) (Å)
- `F_equilibrium`: Minimum F at each T (eV/atom)

### `outfile.ac_equilibrium.dat` (a-c grid)
### `outfile.act_equilibrium.dat` (a-c-T grid)

Plain text file with equilibrium lattice parameters vs temperature:
```
# Temperature (K)    a_eq (Angstrom)    c_eq (Angstrom)    F_min (eV/atom)
1.000000E+00         2.950000E+00       4.650000E+00       -5.123456E+00
2.020202E+01         2.951234E+00       4.651234E+00       -5.122345E+00
...
```

## Theory

### Free Energy Components

The Helmholtz free energy is computed as:

$$
F(x, T) = U(x) + \Delta U_0(x) + F_{ph}(x, T)
$$

where:
- $U(x)$ = Static internal energy from equation of state
- $\Delta U_0(x)$ = Correction to static energy from force constant fitting
- $F_{ph}(x, T)$ = Phonon free energy from quasi-harmonic approximation

### Phonon Free Energy

The phonon free energy per atom is:

$$
F_{ph} = \frac{1}{N_a} \sum_{\mathbf{q},\nu} w_{\mathbf{q}} \left[ \frac{\hbar\omega_{\mathbf{q}\nu}}{2} + k_B T \ln\left(1 - e^{-\hbar\omega_{\mathbf{q}\nu}/k_B T}\right) \right]
$$

where:
- $\omega_{\mathbf{q}\nu}$ = Phonon frequency at wavevector $\mathbf{q}$, branch $\nu$
- $w_{\mathbf{q}}$ = Integration weight
- $N_a$ = Number of atoms

### 4th Order Polynomial Fitting (a-c grid)

For a-c grids, the free energy surface is fitted to a 4th order polynomial:

$$
F(a, c) = \sum_{i+j \leq 4} C_{ij} \tilde{a}^i \tilde{c}^j
$$

where $\tilde{a}$ and $\tilde{c}$ are normalized coordinates for numerical stability.

The equilibrium lattice parameters are found by minimizing this surface.

## Tips and Best Practices

1. **Grid density**: Use at least 5 points per dimension for reliable interpolation

2. **q-mesh convergence**: Test convergence of phonon free energy with q-mesh density

3. **Unstable modes**: If imaginary frequencies occur, those points are flagged (F = 1.234×10⁸) and excluded from fitting

4. **Temperature range**: Default is 1-2000 K with 100 points. Modify source code for different ranges.

5. **Memory**: Large grids with fine q-meshes require significant memory. Use MPI for parallel execution.

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "NOT DONE" error | Check evalmode in infile.evalpoints matches grid type |
| Large F values | Imaginary frequencies present - check structural stability |
| Missing files | Ensure all required input files are present |
| Memory error | Reduce q-mesh density or use more MPI ranks |

## Related Programs

- `fitmultipole`: Fit force constants across parameter grid
- `extract_forceconstants`: Extract force constants from MD trajectories
- `phonon_dispersion_relations`: Compute phonon dispersions
- `thermal_conductivity`: Compute thermal conductivity

## References

1. O. Hellman et al., "Temperature dependent effective potential method for accurate free energy calculations of solids", Phys. Rev. B 87, 104111 (2013)

2. O. Hellman and I. A. Abrikosov, "Temperature-dependent effective third-order interatomic force constants from first principles", Phys. Rev. B 88, 144301 (2013)
