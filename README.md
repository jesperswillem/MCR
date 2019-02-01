# MCR

These scripts use dask and will use all cores available, it is recommended that you run them on a server using a workload manager like slurm.

Setting up a conda environment on server:
>module load miniconda/3
>conda create -c rdkit -n miniconda_rdkit rdkit dask matplotlib scikit-learn
>source activate miniconda_rdkit

Replace the module load command with whateve is applicable to your server. If the server runs conda 2 add python=3 to the conda create command, this code was written and tested for python 3 and might not run in python 2. 

Running a substructure search:
>python search.py CC1CCC1CC=O -d path/to/smiles.smi -o ./output_file.smi

Here it will output molecules that exactly match the query, if you want to do similarity search use:
>python search.py CC1CCC1CC=O -d path/to/smiles.smi -o ./output_file.smi -s 0.95

This will calculate the tanomoto similarity between the topological fingerprints and apply a cut-off of 95% smilarity.
Using a substructure match is however about 6 to 7 times faster.

To view the command line interface (CLI):
>python search.py -h

Running multi-component "reaction" (MCR):
>python MCS.py "C([*:2])C1[C@H](\[*:1])NC(=O)NC=1C" -d path/to/smiles.smi -q "[CX3H1](=O)" "[CH2][CX3]([CH3])(=O)" -m 25

Here the order of the query arguments matter, they have to match wildcards in the scaffold (first argument).

Again for the full CLI:
>python MCR.py -h
