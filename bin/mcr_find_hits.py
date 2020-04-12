#!/usr/bin/env python3
import argparse
import os
import sqlite3
import configparser
import numpy as np

from pathlib import Path
from scipy.spatial import distance
from collections import namedtuple

from rdkit import Chem
from rdkit.Chem import Draw

from MCR import chem_functions
from MCR import helpers
from MCR import core
from MCR.sql import get_chembl_activies


def parse_inputs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--cluster_dir",
                        dest="cluster_dir",
                        default=None,
                        help="Directory containing clusters and centroids...")

    parser.add_argument("-o", "--out_file",
                        dest="out_file",
                        default="cluster_matching.out",
                        help="Outfile to write matching or close chemical matter to.")

    parser.add_argument("-db", "--database",
                        dest="database",
                        default=None,
                        help="Path to SQLite .db file of the ChEMBL database.")

    parser.add_argument("-af", "--activity_file",
                        dest="activity_file",
                        default=None,
                        help="File containing tab-separated activities following the following layout:\n"
                             "target_id*, target_name*, compound_id, canonical_smiles*, standard_type, "
                             "standard_relation, pchembl_value*, standard_units\t (* = mandatory)")

    parser.add_argument("-t", "--target_id",
                        dest="target_id",
                        default=None,
                        help="ChEMBL id of the target for which you want to cluster the known bioactivity data.")

    parser.add_argument('--activity-cutoff',
                        dest="activity_cutoff",
                        default=6.5,
                        help="Cut-off for bioactivity level, for ChEMBL 6.5 seems to be a reasonable number."
                        )

    parser.add_argument('--num_clusters',
                        dest='num_clusters',
                        default=10,
                        help='Number of clusters close to active chemical matter to unpack.')

    options = parser.parse_args()

    return options


def main():
    args = parse_inputs()
    number_cluster_unpacks = 5

    reference_activities = []
    # Read cluster centroids
    centroids = np.loadtxt(args.cluster_dir + 'centroids.txt', delimiter='/t')

    # Query ChEMBL and calculate arrays.
    if args.database and args.target_id:
        print(f'Querying {args.database} for target: {args.target_id}')
        reference_activities.append(get_chembl_activies(args.database, [args.target_id]))

        with open(args.cluster_dir + f'chembl_query.tsv', 'w') as out_file:
            for activity in reference_activities:
                out_file.write('\t'.join(map(str, activity)) + '\n')

    if args.activity_file:
        activity_tuple = namedtuple('activity_tuple',
                                    'target_id, target_name, compound_id, canonical_smiles, standard_type, standard_relation, pchembl_value, standard_units')

        with open(args.cluster_dir + args.activity_file, 'r') as in_file:
            for line in in_file:
                new_activity = line.strip('\n').split('\t')
                new_activity[6] = float(new_activity[6])
                new_activity = activity_tuple(*new_activity)
                reference_activities.append(new_activity)

    reference_activities = [act for act in reference_activities if act.pchembl_value > args.activity_cutoff]

    # Calculate descriptors of activity mols
    if not reference_activities:
        print('Did not retrieve compounds from database or file, exiting.')
        exit(1)

    reference_descriptors = []
    reference_mols = [Chem.MolFromSmiles(act.canonical_smiles) for act in reference_activities]
    reference_mols = [chem_functions.standardize_mol(mol) for mol in reference_mols]
    reference_matrix = chem_functions.generate_fingerprints_iter(reference_mols)

    similar_centroid_matches = {i:(0.0, None, None) for i in range(centroids.shape[0])}
    for reference_index, reference_fp in enumerate(reference_matrix):
        for centroid_index, centroid in enumerate(centroids):
            similarity = 1 - distance.jaccard(reference_fp, centroid)
            if similar_centroid_matches[centroid_index][0] < similarity:
                similar_centroid_matches[centroid_index] = (similarity, reference_index, centroid_index)

    similar_centroid_matches = sorted(similar_centroid_matches.values(), key=lambda x: x[0])
    similar_centroid_matches = [value for value in similar_centroid_matches if value[1] is not None]

    picked_clusteres = [match[2] for match in similar_centroid_matches[:args.num_clusters]]

    # Get cluster files
    cluster_files = []
    for cluster in picked_clusteres:
        cluster_files.append(f"{args.cluster_dir}/MCR_cluster_{cluster}.txt.gz")

    # Process input file and find reactants, and scaffolds.
    config = configparser.ConfigParser()
    config.read(f"{args.cluster_dir}/cluster.inp")

    scaffold = config.get('cluster_info', 'scaffold')
    reactant_files = []
    reacting_groups = []

    for i in range(config.getint('cluster_info', 'num_reactant_files')):
        reactant_files.append(f"{args.cluster_dir}{config.get('cluster_info', f'file{i}')}")
        reacting_groups.append(config.get('cluster_info', f'reacting_group{i}'))

    # Putting matches (couples of similar compound in the reference and clusters) into a named tuple to make it easier
    # to play with during analysis.
    match = namedtuple(
        'match', 'similarity, reference_index, reference_mol, reference_data, cluster, cluster_index, cluster_mol, '
                  'cluster_reactants'
    )
    best_matches = []

    # Load_reactants and reacting groups
    reactants = []
    for path in reactant_files:
        reactant = []
        with open(path, 'r') as f:
            lines = f.readlines()
            compounds = [line.strip('\n').split('\t') for line in lines if 'canonical_smiles' not in line]
            for compound in compounds:
                mol = Chem.MolFromSmiles(compound[0])
                if not mol:
                    continue
                reactant.append(mol)
        reactants.append(reactant)

    # Generate cluster using function in core and write to file.
    for file in cluster_files:

        cluster_index = file.split("/")[-1]
        print(f'Searching {cluster_index}')
        cluster = core.construct_cluster_rxn(file, scaffold, list(zip(reactant_files, reacting_groups)))

        # Generate a feature array
        cluster_mols = [compound[0] for compound in cluster]
        cluster_feature_array = chem_functions.generate_fingerprints_iter(cluster_mols)

        best_similarity = 0.0
        best_reference_index = None
        best_cluster_index = None

        # Searching reference matrix for matches.
        for reference_index, reference_fp in enumerate(reference_matrix):
            for cluster_comp_index, cluster_comp_fp in enumerate(cluster_feature_array):
                similarity = 1 - distance.jaccard(reference_fp, cluster_comp_fp)
                if best_similarity < similarity:
                    best_similarity = similarity
                    best_reference_index = reference_index
                    best_cluster_comp_index = cluster_comp_index

        # 'couple', 'similarity, reference_index, reference_mol, reference_data, cluster_index, cluster_comp_index, cluster_comp_mol'
        reference_mol = reference_mols[best_reference_index]
        reference_data = reference_activities[best_reference_index]
        cluster_comp_mol = cluster_mols[cluster_comp_index]
        match_reactants = []

        for num, reactant_index in enumerate(cluster[cluster_comp_index][2]):
            match_reactants.append(reactants[num][reactant_index])


        new_match = match(
            best_similarity, best_reference_index, reference_mol, reference_data, cluster_index, cluster_comp_index,
            cluster_comp_mol, match_reactants
        )

        best_matches.append(new_match)

    # Finally write out results

    # First make sure we have a directory to write to
    results_dir = f"{args.cluster_dir}/match_result_{str(args.target_id)}/"
    backup = helpers.backup_dir(results_dir)

    Path(results_dir).mkdir()

    # Write best matches to tsv
    matches_tsv = []
    with open(results_dir + 'matches.tsv', 'w') as matches_out:
        for match in best_matches:
            match = [match.similarity, *match.reference_data, match.cluster, match.cluster_index,
                     Chem.MolToSmiles(match.cluster_mol), *[Chem.MolToSmiles(mol) for mol in match.cluster_reactants]]
            match = '\t'.join(map(str, match))
            matches_out.write(match + '\n')

    # Generate rows with reactants, product and reference match.
    match_grid = []
    legends = []
    for match in best_matches:
        match_grid += [*match.cluster_reactants, match.cluster_mol, match.reference_mol]
        legends += ['reactant1', 'reactant2', 'reactant3', str(match.cluster_index),
                   f'{match.reference_index} | '
                   f'target: {match.reference_data.target_name}']

    Draw.MolsToGridImage(mols=match_grid, molsPerRow=5, legends=legends).save(results_dir + 'matches_grid_imagae.png')


if __name__ == '__main__':
    main()
