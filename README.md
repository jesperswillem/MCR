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
conda install -c conda-forge dask dask-ml numpy 
```

Clone this repository
```
git clone https://github.com/jesperswillem/MCR.git
```

And finally run the setup script
```
cd MCR
python setup.py install
```

#### Testing installation

The following command line tools shoud be added to your path after installation
```
construct_clusters.py
mcr_plot.py
run_mcr.py
``` 

To run unittests
```
python -m unittest tests
``` 
<br>
More to follow.

#### Authors
Florian van der Ent<br>
Willem Jespers <br>
Hugo Gutierrez de Teran

