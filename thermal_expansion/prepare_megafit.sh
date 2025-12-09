#!/bin/bash
# Prepare megafit input files for QHA thermal expansion calculation.
# This script runs generate_megafit.py to:
#   1. Run pack_simulation in each a*_c* folder (if needed)
#   2. Extract lattice parameters and DFT energies
#   3. Create infile.simulations and infile.evalpoints
#   4. Print the megafit command to run

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=================================================="
echo "Preparing megafit input files for QHA"
echo "Working directory: $SCRIPT_DIR"
echo "=================================================="
echo

python3 generate_megafit.py

echo
echo "Done! Review the infile.simulations and infile.evalpoints files,"
echo "then run the megafit command printed above."
