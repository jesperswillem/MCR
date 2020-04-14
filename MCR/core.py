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
import copy
from rdkit import Chem
from rdkit.Chem import AllChem


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
    # Hot fix for square brackets with a comma inside.
    for inum, value in enumerate(split):
        if '[' in value and not ']' in value:
            if ']' not in split[inum+1]:
                raise ValueError(f'No closing bracket found in {input}')
            else:
                split[inum] = value + ',' + split[inum+1]
                split.pop(inum + 1)

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
            bag = bag.filter(lambda x, smarts=copy.deepcopy(include_group): x.HasSubstructMatch(smarts))

    if query['exclude_smarts']:
        if type(query['exclude_smarts']) != list:
            query['exclude_smarts'] = [query['exclude_smarts']]
        for exclude_smarts in query['exclude_smarts']:
            exclude_group = Chem.MolFromSmarts(exclude_smarts)
            bag = bag.filter(lambda x, smarts=copy.deepcopy(exclude_group): not x.HasSubstructMatch(smarts))

    return bag


def execute_queries(settings, query_parameters, output_path):
    # TOPOLISH: test and doc
    # Create dask delayed object, make rdkit mol object if possible.
    if settings['reactant_db']:
        input_bag = db.read_text(settings['reactant_db'], blocksize=settings["partition_size"])
        input_bag = input_bag.map(lambda x: Chem.MolFromSmiles(x)).remove(lambda x: x is None)
    else:
        print('Warning no reactant database defined, this only works if all reactants are loaded from sequence or file.')

    # Creating a list of product bags and the compute section.
    reactant_bags = []
    for i, query in enumerate(query_parameters):
        # Query and repartition mol database
        if query["from_file"]:
            bag = db.read_text(query['from_file'])
            bag = bag.map(lambda x: Chem.MolFromSmiles(x)).remove(lambda x: x is None)
        elif query["from_sequence"]:
            bag = db.from_sequence(query["from_sequence"])
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


def execute_mcr(reactants_lists, reacting_groups, scaffolds):
    # TOPOLISH: test and doc
    # TOPOLISH: rdkit reaction version
    for i, reacting_group in enumerate(reacting_groups):
        reacting_group = Chem.MolFromSmarts(reacting_group)
        reactants_lists[i] = [
            [chem_functions.substitute_reactive_group(reactant[0], reacting_group, i + 1)[0], reactant[1]] for
            reactant in reactants_lists[i]]

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

    scaffold = scaffolds[0]

    reactant_product = reactant_product.map(
        lambda x, y=copy.deepcopy(scaffold): chem_functions.join_fragments(x, y))

    mcr_result_bag = reactant_product.map(lambda x: [chem_functions.weld_r_groups(x[0])] + x[1:])

    return mcr_result_bag, expected_total


def execute_mcr_rxn(reactants_lists, reacting_groups, scaffolds):
    # TOPOLISH: test and doc

    # predict number of compound to select a sample later on.
    expected_total = 1
    for k in reactants_lists:
        num_reactants = len(k)
        if num_reactants == 0:
            raise Exception('Empty reactants file')
        expected_total *= len(k)

    #create dask bags
    reactants_lists = [db.from_sequence(k, partition_size=int(10)) for k in reactants_lists]

    # Make a single product bag of the form [['smiles', reactant01_index, reactant02_index], [...], ...]
    reactant_product = reactants_lists.pop(0)
    for p in reactants_lists:
        reactant_product = reactant_product.product(p)
    reactant_product = reactant_product.map(helpers.flatten_tuple)

    # The product come in as [[a, id_a], [b, id_b]] my_map is used to make [[a, b], [id_a, id_b]]
    # So later on we can get [product_a_b, [id_a, id_b]]
    def my_map(x):
        a = []
        b = []
        for y in x:
            a.append(y[0])
            b.append(y[1])
        return a, b

    # Join reactants by dots to produce a rxn reactant later.
    reactants_groups = '.'.join(reacting_groups)

    # Loop over scaffolds and perform reactions
    for scaffold in scaffolds:
        rxn_reaction = f"{reactants_groups}>>{scaffold}"
        reactant_product = reactant_product.map(my_map)
        mcr_result_bag = reactant_product.map_partitions(lambda x, rxn=rxn_reaction: chem_functions.apply_rxn(x, rxn))
    mcr_result_bag = mcr_result_bag.flatten()

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
    # TODO: add test/doc culling empty partitions twice now, also happens in molbag to featurearray.
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


def construct_cluster_rxn(cluster_file, scaffold, reactant_data):
    """recreates a cluster from the original data if that cluster has been created using rxn.

    Parameters
    ----------
    cluster_file: str
    scaffold: str
    reactant_data: list(list)

    Returns
    -------
    cluster: list(tuple)
    """
    # Load_reactants and reacting groups
    reacting_groups = []
    reactants = []
    for reacting_group in reactant_data:

        reacting_groups.append(reacting_group[1])
        reactant = []
        with open(reacting_group[0], 'r') as f:
            lines = f.readlines()
            compounds = [line.strip('\n').split('\t') for line in lines if 'canonical_smiles' not in line]
            for compound in compounds:
                mol = Chem.MolFromSmiles(compound[0])
                if not mol:
                    continue
                reactant.append(mol)
        reactants.append(reactant)

    # Unzipping cluster files and constructing list from lines
    with gzip.open(cluster_file, 'rt') as f:
        lines = f.readlines()

    mols_to_rebuild = [line.strip('\n').split('\t') for line in lines]

    # Create rdkit reaction from smart
    rxn_reaction = f"{'.'.join(reacting_groups)}>>{scaffold}"
    reaction = AllChem.ReactionFromSmarts(rxn_reaction)

    # Perform reaction and cleaup products.
    cluster = []
    for mol_to_rebuild in mols_to_rebuild:
        # Note to self: since I stored python lists as string I use eval to recreate lists. In other software this might
        # be a major security risk.
        building_block_indices, compound_index = mol_to_rebuild
        building_block_indices = eval(building_block_indices)
        building_blocks = [reactants[index][block] for index, block in enumerate(building_block_indices)]
        unfiltered_products = reaction.RunReactants(building_blocks)
        found_smiles = set()
        for product in unfiltered_products:
            product = product[0]
            smiles = Chem.MolToSmiles(product)
            if not smiles in found_smiles:
                # Chem.GetSSSR(product)
                Chem.SanitizeMol(product)
                cluster.append((product, compound_index, building_block_indices))
                found_smiles.add(smiles)

    return cluster


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
