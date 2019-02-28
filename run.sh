#!/bin/bash
#sbatch -n 8  -t 6:00:00 -J jobname script.sh

module load miniconda/3
source activate miniconda_rdkit

#python search.py CC1CCC1CC=O -d ~/databases/emolecules/emol_100000.smi -o ./output_file.smi -s 0.95
#python MCS.py "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C" -d ~/databases/emolecules/emol_100000.smi -q "[CX3H1](=O)" "[CH2][CX3]([CH3])(=O)" -m 25
#python cluster_compounds.py -o output_clustering.png -db /home/vonderent/databases/ChEMBL/chembl_24_sqlite/chembl_24.db -t CHEMBL1833 CHEMBL255 -pca -i /home/vonderent/1.A2b_ago/0.lig_gen_and_docking/1.clustering/generated_outfile.smi
python QSAR.py -db /home/vonderent/databases/ChEMBL/chembl_24_sqlite/chembl_24.db -t CHEMBL256 -i /home/vonderent/1.A2b_ago/0.lig_gen_and_docking/1.clustering/generated_outfile.smi -o compound_predictions.smi
