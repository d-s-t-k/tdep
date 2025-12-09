#!/usr/bin/env python3
import os
from ase.io import read
from ase.calculators.espresso import Espresso, EspressoProfile

# --- User settings ---
infile = "infile.ucposcar"   # VASP-format input
run_dir = "./ev"             # QE run directory

# QE setup
pseudopotentials = {
    "Li": "li_pbesol_v1.4.uspp.F.UPF",
    "Nb": "Nb.pbesol-spn-kjpaw_psl.0.3.0.UPF",
    "O":  "O.pbesol-n-kjpaw_psl.0.1.UPF",
}
input_data = {
    "control": {"calculation":"scf","tprnfor":True,"tstress":True,
                "verbosity":"medium","disk_io":"none"},
    "system":  {"ecutwfc":60,"ecutrho":400,"occupations":"fixed","nosym":False},
    "electrons":{"conv_thr":1.0e-12},
}
profile = EspressoProfile(
    command="pw.x",
    pseudo_dir="/home/qhossain/SSSP_1.3.0_PBEsol_efficiency"
)

# --- Prepare run directory ---
os.makedirs(run_dir, exist_ok=True)

# --- Read structure and run QE ---
atoms = read(infile, format="vasp")
atoms.pbc = [True, True, True]

calc = Espresso(
    profile=profile,
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    directory=run_dir,        # QE runs here; writes pw.in / pw.out
    kpts=(6, 6, 6),
    koffset=(0, 0, 0),
)
atoms.calc = calc
atoms.get_potential_energy()  # triggers QE, writes pw.out

print(f"Done. Check {run_dir}/pw.out")
