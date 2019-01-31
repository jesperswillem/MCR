#!/bin/bash
#sbatch -n 8  -t 6:00:00 -J jobname script.sh

module load miniconda/3
source activate miniconda_rdkit
python search.py CC1CCC1CC=O -d /home/vonderent/databases/emolecules/emol_full.smi -o ./output_file.smi -s 0.8
