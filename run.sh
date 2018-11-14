#!/bin/bash
#sbatch -n 8  -t 6:00:00 -J jobname script.sh

module load miniconda/3
source activate miniconda_rdkit
python main.py run
