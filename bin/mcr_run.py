#!/usr/bin/env python3
""" Command line script of that comes with this library, takes an input file and outputs a query, reactants resulting of a
MCR or both.
"""
import argparse
import configparser
import copy
import sys
from pathlib import Path

import dask
from dask.diagnostics import ProgressBar
from rdkit import RDLogger, rdBase, Chem

from MCR import helpers, chem_functions, core

def settings_print(x):
    """Prints dict x in a readable format.

    Parameters
    ----------
    x: dict

    Returns
    -------
    """
    for k, v in x.items():
        padding = 55 - len(str(k)) - len(str(v))
        print(f"{k}{padding*' '}{v}")


def startup_report(settings, query_parameters, MCR_parameters):
    """Takes in the processed data from the input file and prints the settings that will be used.

    Parameters
    ----------
    settings: dict
    query_parameters: list[dict]
    MCR_parameters: dict

    Returns
    -------
    """
    print("------------------Used-settings------------------\n")
    settings_print(settings)

    for i, parameters in enumerate(query_parameters):
        print(f"\nParameters for reactant {i+1}:")
        settings_print(parameters)

    if MCR_parameters['perform_mcr']:
        print("\nSettings for multi-component reaction:")

        settings_print(MCR_parameters)

    print("\n------------------------------------------------\n")


def parse_input():
    """Function that reads in the comment line arguments, if a input file is given it will parse the input file and
    return MCR/query parameters.

    Parameters
    ----------

    Returns
    -------
    settings: dict
    query_parameters: list(dict, ..., m, n)
    MCR_parameters: dict
    """

    # Parse the command line input to get the input file
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file",
                        help="Path to input file giving instructions for queries and/or MCR reaction.")
    options = parser.parse_args()

    # Make sure an input file was parsed
    if not options.input_file:
        print(">>> An input file must be specified.")
        parser.print_help()
        sys.exit(1)

    # Process config file and construct dictionaries used to construct pipeline
    config = configparser.ConfigParser()
    config.read(options.input_file)
    with open(options.input_file) as f:
        config.read(f)

    # Create settings dictionary with standard settings.
    settings = {
        "num_processes": 1,
        "partition_size": 1e6,
        "reactant_db": None,
        "log_file": "MCR.log",
        "output_folder": None
    }

    if config.has_section("general"):
        for option in config["general"]:
            if option not in settings:
                raise KeyError(f"Invalid option {option} in general section.")
            settings[option] = core.parse_input_field(config["general"][option])


    # Create new_reactant dictionary with default settings
    new_reactant_defaults = {
        "from_file": None,
        "from_sequence": None,
        "reacting_group": None,
        "include_smarts": None,
        "exclude_smarts": None,
        "max_heavy_ratoms": 10,
        "wash_molecules": True,
        "keep_largest_fragment": False
    }

    reactants = []
    for i in range(1, 10):
        j = "reactant" + str(i).zfill(2)
        new_reactant = copy.deepcopy(new_reactant_defaults)
        if config.has_section(j):
            for option in config[j]:
                if not option in new_reactant:
                    raise KeyError(f"Invalid option {option} in {j} section.")
                new_reactant[option] = core.parse_input_field(config[j][option])
            try:
                if new_reactant["reacting_group"] != None:
                    reactants.append(new_reactant)
                else:
                    raise RuntimeError
            except:
                raise ValueError(f"No smarts defined in {j}")
        else:
            break

    # Create new_reactant dictionary with default settings
    MCR_parameters = {
        "perform_mcr": False,
        "scaffold": None,
        "subset_method": None,
        "subset_size": 1e4,
        "training_sample": 1e5,
        "max_heavy_ratoms": 35,
        "remove_pains": False,
        "remove_non_lipinski": False,
    }

    if config.has_section("MCR"):
        for option in config["MCR"]:
            if not option in MCR_parameters:
                raise KeyError(f"Invalid option {option} in MCR section.")
            MCR_parameters[option] = core.parse_input_field(config["MCR"][option])
        if MCR_parameters["scaffold"]:
            if type(MCR_parameters["scaffold"]) == str:
                MCR_parameters["scaffold"] = [MCR_parameters["scaffold"]]
            MCR_parameters["perform_mcr"] = True

    return settings, reactants, MCR_parameters


def write_config(path, settings, query_parameters, MCR_parameters):
    """Writes out settings dict to .inp file at given path.

    Parameters
    ----------
    path: str
    settings: dict
    query_parameters: dict
    MCR_parameters: dict

    Returns
    -------
    """

    config = configparser.ConfigParser()

    config.add_section('general')
    for key, value in settings.items():
        config.set("general", key, str(value))

    for i,query in enumerate(query_parameters):
        section = f'reactant{str(i).zfill(2)}'
        config.add_section(section)
        for key, value in query.items():
            config.set(section, key, str(value))

    config.add_section('MCR')
    for key, value in MCR_parameters.items():
        config.set("MCR", key, str(value))

    with open(path, 'w') as f:
        config.write(f)


def write_products(bag, path, file_name):
    """Write out a bag of mols with indices of reactants to file
    """
    bag = bag.map(lambda product: Chem.MolToSmiles(product[0]) + '\t' + '\t'.join(map(str, product[1])) + '\n')
    # bag.to_textfiles(path + file_name) # Didn't use this because it splits the file (it's faster though)
    with ProgressBar():
        with open(path + file_name, 'w') as f:
            for product in bag.compute():
                f.write(product)


def main():
    """Main function for running from the command line with an input file.
    """

    RDLogger.DisableLog('rdApp.info')
    rdBase.DisableLog('rdApp.error')

    # Read input file into dictionaries used later.
    settings, query_parameters, MCR_parameters = parse_input()
    startup_report(settings, query_parameters, MCR_parameters)

    # Make sure we have a output directory to write to.
    if not settings['output_folder']:
        print('No output directory defined for MCR!')
        sys.exit(1)
    output_path = settings['output_folder']
    if output_path[-1] != "/":
        output_path += "/"

    # Backup existing folder with the same name if needed.
    backup = helpers.backup_dir(output_path[:-1])
    if backup:
        print(f"Backed up existing {output_path[:-1]} folder to {backup}\n")
    Path(output_path).mkdir(parents=True, exist_ok=True)
    output_path = str(Path(output_path).absolute())

    if output_path[-1] != "/":
        output_path += "/"

    # Write out used input file.
    write_config(output_path + "used_input_file.inp", settings, query_parameters, MCR_parameters)

    # Execute and write out queries
    reactants_lists = core.execute_queries(settings, query_parameters, output_path)

    # Check if MCR is requested
    if MCR_parameters['perform_mcr']:
        reacting_groups = [query['reacting_group'] for query in query_parameters]
        # Old manual reaction function
        # mcr_result, expected_total = core.execute_mcr(reactants_lists, reacting_groups, MCR_parameters['scaffold'])

        # rdkit implementation
        mcr_result, expected_total = core.execute_mcr_rxn(reactants_lists, reacting_groups, MCR_parameters['scaffold'])

        if type(MCR_parameters['max_heavy_ratoms']) == int:
            mcr_result = mcr_result.filter(
                lambda x, y=copy.deepcopy(MCR_parameters['max_heavy_ratoms']): x[0].GetNumHeavyAtoms() <= y)
        else:
            raise ValueError

        if MCR_parameters["remove_non_lipinski"]:
            mcr_result = mcr_result.filter(lambda x: chem_functions.filter_lipinski(x[0]))

        if MCR_parameters["remove_pains"]:
            mcr_result = mcr_result.filter(lambda x: chem_functions.filter_pains(x[0]))

        # Perform clustering if requested
        if MCR_parameters['subset_method'] == "k-means":
            print('Subsampling using k-means clustering.')
            # Create training set to find for KMeans
            sample_prob = MCR_parameters["training_sample"]/expected_total
            centroids, training_sample = core.make_kmeans_clusters(mcr_result, sample_prob)

            # Assign clusters training sample, kind of useless the centroids are not real points.
            # train_desc_bag = training_sample.map(lambda x: [chem_functions.generate_fingerprints(x), Chem.MolToSmiles(x)])
            # train_assigned_node = train_desc_bag.map(lambda x, centers=centroids: [helpers.match_node(x[0], centers)] + [x[1:]])

            # Assign clusters
            desc_bag = mcr_result.map(lambda x: [chem_functions.generate_fingerprints(x[0])] + list(x[1:]))
            prediction = desc_bag.map(lambda x, centers=centroids: [helpers.closest_node(x[0], centers)] + [x[1:]])

            # Compute by creating delayed cluster writers for each partition
            delayed_partitions = prediction.to_delayed(optimize_graph=True)

            writers = []
            for partition_index, partition in enumerate(delayed_partitions):
                writers.append(dask.delayed(lambda x, y=copy.deepcopy(partition_index): helpers.cluster_writer(
                    output_path + "MCR", y, x))(partition))

            print('\nAssigning labels and writing full set to disk...')
            with ProgressBar():
                dask.compute(*writers)

            # Combine the cluster files from different partitions into one file per cluster and enumerate the compounds
            # in each cluster
            print('\nCombining partitions into single file per cluster...')
            helpers.combine_files(output_path + "MCR")

            # Write out a inp file used to reconstruct the clusters from .txt.gz files.
            cluster_config = configparser.ConfigParser()
            cluster_config.add_section('cluster_info')

            cluster_config.set('cluster_info', 'num_reactant_files', str(len(reactants_lists)))
            for partition_index, query in enumerate(query_parameters):
                cluster_config.set('cluster_info', f'file{partition_index}', f"./reactant{str(partition_index).zfill(2)}.smi")
                cluster_config.set('cluster_info', f'reacting_group{partition_index}', query['reacting_group'])

            cluster_config.set('cluster_info', 'scaffold', ', '.join(MCR_parameters['scaffold']))

            helpers.store_np_matrix(centroids, f'{output_path}centroids.txt')

            with open(f'{output_path}cluster.inp', 'w') as f:
                cluster_config.write(f)
        elif MCR_parameters['subset_method'].lower() in ['false', 'no', 'not', '0']:
            write_products(mcr_result, output_path, 'mcr_result.smi')
        else:
            print(f"Error: unkown subsampling method: {MCR_parameters['subset_method']}, recognised options are: k-means or False")
            print("Writing out products without subsampling...")
            write_products(mcr_result, output_path, 'mcr_result.smi')


if __name__ == "__main__":
    main()
