#!/usr/bin/env python3
import argparse
import sys
import configparser
from pathlib import Path
from MCR import core
from rdkit import rdBase, RDLogger

def read_inputs():
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
                        help="Path to input file giving instructions on the used reactant files and scaffold.")

    parser.add_argument("clusters",
                        nargs='+',
                        type=int,
                        help="MCR cluster index, at least one needs to be provided.")

    options = parser.parse_args()

    # Make sure an input file was parsed
    if not options.input_file:
        print(">>> An input file must be specified.")
        parser.print_help()
        sys.exit(1)

    # Get directory to read form/write to
    directory = Path(options.input_file).parts[:-1]
    directory = '/'.join(directory)

    # Process input file and find reactants, and scaffolds.
    config = configparser.ConfigParser()
    config.read(options.input_file)

    reactant_files = []
    reacting_groups = []
    for i in range(config.getint('cluster_info', 'num_reactant_files')):
        reactant_files.append(f"{directory}/{config.get('cluster_info', f'file{i}')}")
        reacting_groups.append(config.get('cluster_info', f'reacting_group{i}'))

    scaffold = config.get('cluster_info', 'scaffold')

    # Get cluster files
    cluster_files = []
    for cluster in options.clusters:
        cluster_files.append(f"{directory}/MCR_cluster_{cluster}.txt.gz")

    return cluster_files, scaffold, list(zip(reactant_files, reacting_groups)), directory


def main():
    # Suppress warnings for headers.
    RDLogger.DisableLog('rdApp.info')
    rdBase.DisableLog('rdApp.error')

    cluster_files, scaffold, reactants, directory = read_inputs()

    # Generate cluster using function in core and write to file.
    for file in cluster_files:
        cluster = core.construct_cluster(file, scaffold, reactants)
        outfile = f"{directory}/constructed_{file.split('/')[-1].split('.')[0]}.smi"
        with open(outfile, 'w') as f:
            f.writelines('\n'.join(cluster))


if __name__ == '__main__':
    main()
