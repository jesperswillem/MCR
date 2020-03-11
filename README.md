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


#### Testing installation

To run unittests
```
cd MCR/tests
python -m unittest
``` 
<br>
More to follow.


#### Authors
Willem Jespers <br>
Florian van der Ent<br>
Hugo Gutierrez de Teran
