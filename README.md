Set of tools and a library for the generation of virtual screening libraries from multi-component reactions.

#### Requirements

Currently requires:
- Dask
- Dask-ml
- rdkit
- numpy


#### Installation (Linux)

Create a new conda env
```
conda create -c rdkit -n MCR_env rdkit
```

Install required packages
```
source activate MCR_env
conda install dask dask-ml numpy 
```

Clone this repository or the devel branch
```
# master
git clone https://github.com/jesperswillem/MCR.git
# devel
git clone -b devel --single-branch https://github.com/jesperswillem/MCR.git 
```

And then install by running the setup script
```
cd MCR
python setup.py install
```
The following command line tools should now have been added to your path
```
construct_clusters.py
mcr_plot.py
run_mcr.py
``` 
If these do not show up in your path try reloading your conda environment after installation.

#### Testing installation

To run unittests
```
cd MCR/tests
python -m unittest
``` 
<br>

#### Tutorial

Two example input files are available in the example folder, in this tutorial we'll have a look at the 
pyrazine-2(1H)-one multi-component reaction. This multi-component reaction is discussed in more detail 
[here](https://pubs.acs.org/doi/abs/10.1021/jo4003163). In this example we reproduce the compounds of the SAR performed 
[here](https://www.future-science.com/doi/abs/10.4155/fmc.15.69?rfr_dat=cr_pub%3Dpubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&journalCode=fmc). 
Before we start generating compounds have a look at the input file (examples/pyrazin-2(1H)-one.inp). A MCR input file 
consists of 3 types of sections. The first is the general section:
```
[general]
partition_size = 2.5MB
output_folder = ./single_step_ugi_from_sequence
``` 
In this case we only specify the output folder here. Next are the reactants section let's have a look at [reactant01]:

```
[reactant01]
reacting_group = [NX3H2:1]
from_sequence = Nc1ccccc1, NCCCCC, NCCN(C)C
```

Reactant one can is an amine we selected 3 amines from the mentioned SAR paper and supply them form sequence. It will 
often interesting to query a database from reagents but when investigating a new reaction it is recommended to first
reproduce existing compounds to validate your input file. 

An important (and often tricky setting) is the reacting group
of each reactant this should be a [SMARTS string](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) with 
at least one numbered atom (in this case the nitrogen which is number one) You should be very specific about the atoms 
in your reactant group. Here we specify a nitrogen with 3 bonds (X3) and 2 hydrogens.

Lastly we will look at the [MCR] section. Here specify the scaffold of our reaction and how we want to process our 
products:

```
[MCR]
scaffold = [C:3]1(=C([N:1](C([C:2](=N1))=O))C([N]([H])[*:4])=O)
subset_method = False
max_heavy_ratoms = 40 !default is 30
``` 

In the scaffold we define where our residues will and up using numbered atoms corresponding to the numbered atoms in our 
reagents. We can also specify filters, here we set the default maximun of 30 heavy atoms to 40. And lastly there is a 
subset_method parameter, This can be either False or k-means. If set to false the products will be written to a single 
.smi file, this is not recommended for input files that generate a large number of products.

To run our MCR we type the following command:
```
mcr_run.py examples/pyrazin-2(1H)-one.inp 
```

This will generate a directory with a copy of the input file and .smi files for the used reagents. In this case where we 
don't cluster our products there will also be a mcr_result.smi file.

#### Authors
Willem Jespers <br>
Florian van der Ent<br>
Hugo Gutierrez de Teran
