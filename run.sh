#!/bin/bash
#sbatch -n 8  -t 6:00:00 -J jobname script.sh

module load miniconda/3
source activate miniconda_rdkit
python MCS.py "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C" -d ~/databases/emolecules/emol_100000.smi -q "[CX3H1](=O)" "[CH2][CX3]([CH3])(=O)" -m 25

