"""Script containing core functionality functions of uppMCR script.
"""

__version__ = "0.2.0"

# Standard library imports
import gzip
import re
import copy
import warnings

from dask_ml import cluster

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Import other parts of this module.
from MCR import chem_functions, helpers

# Imports from other modules
import dask
import dask.bag as db
import dask.array as da
from dask.diagnostics import ProgressBar
import numpy as np
from rdkit import Chem


def parse_input_field(input):
    """Prints dict x in a readable format.

    Parameters
    ----------
    input: str

    Returns
    -------
    output: str or list[str]
    """
    # remove comments
    input = re.split(r'!\s*(?=[^)]*(?:\(|$))', input)[0]

    # split on commas outside of brackets
    split = re.split(r',\s*(?=[^)]*(?:\(|$))', input)
    if len(split) > 1:
        result = []
        for entry in split:
            result.append(re.split(r'!\s*(?=[^)]*(?:\(|$))', entry)[0].strip())
    else:
        result = split[0].strip()
        if type(result) == str:
            try:
                result = int(result)
            except:
                pass
    return result


def cull_empty_partitions(bag):
    """When bags are created by filtering or grouping from a different bag,
    it retains the original bag's partition count, even if a lot of the
    partitions become empty.
    Those extra partitions add overhead, so it's nice to discard them.
    This function drops the empty partitions.

    Parameters
    ----------
    bag: dask.bag

    Returns
    -------
    partitions: dask.bag
    """
    bag = bag.persist()

    def get_len(partition):
        # If the bag is the result of bag.filter(),
        # then each partition is actually a 'filter' object,
        # which has no __len__.
        # In that case, we must convert it to a list first.
        if hasattr(partition, '__len__'):
            return len(partition)
        return len(list(partition))

    partition_lengths = bag.map_partitions(get_len).compute()

    # Convert bag partitions into a list of 'delayed' objects
    lengths_and_partitions = zip(partition_lengths, bag.to_delayed())

    # Drop the ones with empty partitions
    partitions = (p for l, p in lengths_and_partitions if l > 0)

    # return list of delayed objects
    return db.from_delayed(partitions)


def construct_query(query, bag):
    """Takes a dict called query which contains the following parameters:
    reacting_group
    include_smarts
    exclude_smarts
    max_heavy_ratoms
    keep_largest_frag
    wash_molecules

    Then queries the database file(s) to a lazy result that can be used to write to file efficiently or used as a
    component in an MCR graph.

    Parameters
    ----------
    query: dict
    database: str or list[str]

    Returns
    -------
    bag: dask.bag
    """

    # Make rdkit mol object if possible. THis happens in the main module now, let see if that works.
    # bag = bag.map(lambda x:Chem.MolFromSmiles(x )).remove(lambda x: x is None)

    if query['wash_molecules']:
        # Washing protocol:
        # step 1: remove single atoms + some bigger ions
        bag = bag.map(chem_functions.remove_salts_mol)
        # Step 2: Decharge
        bag = bag.map(chem_functions.decharge_mol)
        # Step 3: Standardize
        bag = bag.map(chem_functions.standardize_mol)

    # Step 4?: either keep the largest fragment or remove all mols with more than one fragment.
    if query['keep_largest_fragment']:
        bag = bag.map(chem_functions.get_largest_fragment_mol)
    else:
        bag = bag.filter(lambda x: not (len(chem_functions.split_mol(x)) > 1))

    # Filter for max number of heavy atoms
    if type(query['max_heavy_ratoms']) == int:
        bag = bag.filter(lambda x: x.GetNumHeavyAtoms() <= query['max_heavy_ratoms'])
    else:
        raise ValueError

    # Apply smarts queries.
    if query['reacting_group']:
        reacting_group = Chem.MolFromSmarts(query['reacting_group'])
        bag = bag.filter(lambda x: x.HasSubstructMatch(reacting_group))

    if query['include_smarts']:
        if type(query['include_smarts']) != list:
            query['include_smarts'] = [query['include_smarts']]
        for include_smarts in query['include_smarts']:
            include_group = Chem.MolFromSmarts(include_smarts)
            bag = bag.filter(lambda x: x.HasSubstructMatch(include_group))

    if query['exclude_smarts']:
        if type(query['exclude_smarts']) != list:
            query['exclude_smarts'] = [query['exclude_smarts']]
        for exclude_smarts in query['exclude_smarts']:
            exclude_group = Chem.MolFromSmarts(exclude_smarts)
            bag = bag.filter(lambda x: not x.HasSubstructMatch(exclude_group))

    return bag


def execute_queries(settings, query_parameters, output_path):
    # TOPOLISH: test and doc
    # Create dask delayed object, make rdkit mol object if possible.
    input_bag = db.read_text(settings['reactant_db'], blocksize=settings["partition_size"])
    input_bag = input_bag.map(lambda x: Chem.MolFromSmiles(x)).remove(lambda x: x is None)

    # Creating a list of product bags and the compute section.
    reactant_bags = []
    for i, query in enumerate(query_parameters):
        # Query and repartition mol database
        if query["from_file"]:
            bag = db.read_text(query['from_file'])
            bag = bag.map(lambda x: Chem.MolFromSmiles(x)).remove(lambda x: x is None)
        else:
            bag = construct_query(query, input_bag)
        reactant_bags.append(bag)

    # Compute and enumerate.
    print('Querying database and/or reading input files...')
    with ProgressBar():
        reactants_lists = dask.compute(*reactant_bags)
        # reactants_lists = [common.zfill_enum(i, 6) for i in reactants_lists]
        reactants_lists = [helpers.enum(i) for i in reactants_lists]

    # Write reactants to file.
    print('Writing reactants to file...')
    written = []
    for i, reactants_list in enumerate(reactants_lists):
        output_file = output_path + f"reactant{str(i).zfill(2)}.smi"
        written.append(output_file)
        with open(output_file, 'w') as f:
            f.write("canonical_smiles\tindex")
            for line in reactants_list:
                f.write(f"\n{Chem.MolToSmiles(line[0])}\t{line[1]}")

    return reactants_lists


def execute_mcr(reactants_lists, query_parameters, mcr_parameters):
    # TOPOLISH: test and doc
    for i, query in enumerate(query_parameters):
        if query['reacting_group']:
            reacting_group = Chem.MolFromSmarts(query["reacting_group"])
            reactants_lists[i] = [
                [chem_functions.substitute_reactive_group(reactant[0], reacting_group, i + 1)[0], reactant[1]] for
                reactant in reactants_lists[i]]
        else:
            raise IOError(
                f'Reacting group not defined for reactant {i}, This is required if you want to perform a MCR.')

    expected_total = 1
    for k in reactants_lists:
        expected_total *= len(k)

    reactants_lists = [db.from_sequence(k, partition_size=int(10)) for k in reactants_lists]

    # Make a single product bag of the form [['smiles', reactant01_index, reactant02_index], [...], ...]
    reactant_product = reactants_lists.pop(0)
    for p in reactants_lists:
        reactant_product = reactant_product.product(p)
    reactant_product = reactant_product.map(helpers.flatten_tuple)

    # Make smiles
    reactant_product = reactant_product.map(
        lambda x: [[Chem.MolToSmiles(i[0]), i[1]] for i in x])

    scaffold = mcr_parameters['scaffold']

    reactant_product = reactant_product.map(
        lambda x, y=copy.deepcopy(scaffold): chem_functions.join_fragments(x, y))

    mcr_result_bag = reactant_product.map(lambda x: [chem_functions.weld_r_groups(x[0])] + x[1:])

    # TODO: max heavy atoms.
    # TODO: Lipinski
    # TODO: Lipinski subcomponents
    # TODO: pains filters.

    return mcr_result_bag, expected_total


def molbag_to_featurearray(bag):
    """"Takes a bag of rdkit mol objs and computes an array of morgan fingerprint descriptors of the shape (n, 512).

    Parameters
    ----------
    bag: dask.bag

    Returns
    -------
    features: dask.array
    """
    bag = cull_empty_partitions(bag)
    desc_bag = bag.map_partitions(chem_functions.generate_fingerprints_iter)
    stacked = [da.from_delayed(x, shape=(512, np.NaN), dtype=np.int) for x in
               desc_bag.map_partitions(lambda x: np.stack(x)).to_delayed()]
    [x.compute_chunk_sizes() for x in stacked]
    features = da.concatenate(stacked, axis=0)
    return features


def make_kmeans_clusters(mol_bag, sample_prob):
    # TODO: add test/doc
    training_sample = mol_bag.random_sample(sample_prob).map(lambda x: x[0])
    # Culling empty partitions because they break the generate_fingerprint_iter function.
    training_sample = cull_empty_partitions(training_sample)
    features = molbag_to_featurearray(training_sample).persist()

    # Fit KMeans model and store centroids
    print('\nClustering MCR product training set...')
    model = cluster.KMeans(n_clusters=15, n_jobs=6)
    model.fit(features)
    centroids = model.cluster_centers_

    return centroids, training_sample


def construct_cluster(cluster_file, scaffold, reactants):

    reactant_dicts = []
    for i, j in enumerate(reactants):
        new_dict = {}
        reacting_group = Chem.MolFromSmarts(j[1])
        with open(j[0], 'r') as f:
            lines = f.readlines()
            compounds = [line.strip('\n').split('\t') for line in lines]
            for compound in compounds:
                mol = Chem.MolFromSmiles(compound[0])
                if not mol:
                    continue
                reactant = chem_functions.substitute_reactive_group(mol, reacting_group, i + 1)[0]
                new_dict[compound[1]] = Chem.MolToSmiles(reactant)
        reactant_dicts.append(new_dict)

    with gzip.open(cluster_file, 'rt') as f:
        lines = f.readlines()

    cluster = []
    compounds = [line.strip('\n').split('\t') for line in lines]
    for comp in compounds:
        index = comp.pop()
        to_weld = scaffold
        for i, reactant in enumerate(comp):
            to_weld = f'{to_weld}.{reactant_dicts[i][reactant]}'
        to_weld = Chem.MolFromSmiles(to_weld)
        welded = chem_functions.weld_r_groups(to_weld)
        cluster.append(f'{Chem.MolToSmiles(welded)}\t{index}')

    return cluster
