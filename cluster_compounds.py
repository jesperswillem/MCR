import argparse, sys

import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, SDWriter
from sklearn.manifold import TSNE, MDS, SpectralEmbedding
from sklearn.decomposition import PCA
from collections import namedtuple, Counter

from shared import get_chembl_activies


def cluster_reference_data(targets, mol_activities, methods, output_png, inp_smiles):
    method_dict = {
        "PCA": PCA(n_components=2),
        "t-SNE": TSNE(n_components=2, perplexity=5, learning_rate=40.0, verbose=1),
    }

    plot_dimensions = [7, 7]
    plot_length = plot_dimensions[1] * len(methods)
    plot_width = plot_dimensions[0] * len(targets)
    (fig, subplots) = plt.subplots(len(methods), len(targets), figsize=(plot_width, plot_length))

    # If supplied reference molecules are loaded.
    gen_tuple = namedtuple('gen_tuple', 'target_id, canonical_smiles')
    if inp_smiles:
        with open(inp_smiles) as infile:
            gen_smiles = []
            for line in infile:
                gen_smiles.append(gen_tuple('gen', line.split(' ')[0]))
        mol_activities += gen_smiles
        print("\nRead {} reference molecules".format(len(gen_smiles)))

    # Calculate fingerprints for compounds
    mols = [Chem.MolFromSmiles(x.canonical_smiles) for x in mol_activities]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 4) for m in mols]
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)

    # iterate over chosen methods of clustering to create a row of figures for each
    for k, clus_method in enumerate(methods):
        print('\nPerforming {} clustering...'.format(clus_method))
        X = method_dict[clus_method].fit_transform(np_fps)

        for i, target in enumerate(targets):
            print('Plotting {} activity...'.format(target))
            # get the name for the target being plotted
            for mol in mol_activities:
                if mol.target_id == target:
                    target_name = mol.target_name
                    break
            # find compounds with pChEMBL value for this target so we can colour them using this value
            colors = []
            reds = []
            uncolored = []
            colored = []
            for j, mol in enumerate(mol_activities):
                if mol.target_id == 'gen':
                    reds.append(j)
                elif mol.target_id == target:
                    colors.append(mol.pchembl_value)
                    colored.append(j)
                else:
                    uncolored.append(j)
            colors = np.array(colors)

            # plot compound colouring them based on pchembl value for this target or grey if not available.
            if len(methods) == 1:
                subplots[i].scatter(X[uncolored, 0], X[uncolored, 1], c="grey", marker='.')
                subplots[i].scatter(X[reds, 0], X[reds, 1], c="red", marker='.', alpha=0.5)
                subplots[i].scatter(X[colored, 0], X[colored, 1], c=colors, marker='.', alpha=0.5, cmap='viridis')
                subplots[i].set_title(clus_method + "_" + target_name)
            else:
                subplots[k][i].scatter(X[uncolored, 0], X[uncolored, 1], c="grey", marker='.')
                subplots[k][i].scatter(X[reds, 0], X[reds, 1], c="red", marker='.', alpha=0.5)
                subplots[k][i].scatter(X[colored, 0], X[colored, 1], c=colors, marker='.', alpha=0.5, cmap='viridis')
                subplots[k][i].set_title(clus_method + "_" + target_name)
    plt.savefig(output_png)


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--out_file",
                        dest="out_file",
                        default="clustered.png",
                        help="outfile to write image of clustered compounds to. Make sure to use the .png extension.")

    parser.add_argument("-i", "--input_smiles",
                        dest = "input_smiles",
                        default = None,
                        help = "Input smiles to be plotted..")

    parser.add_argument("-db", "--database",
                        dest="database",
                        default=None,
                        help="path to SQLite .db file of the ChEMBL database.")

    parser.add_argument("-t", "--target_ids",
                        dest="target_ids",
                        default=None,
                        nargs='*',
                        help="ChEMBL ids of the targets for which you want to cluster the known bioactivity data. Requires at least 2 arguments.")

    parser.add_argument("-pca",
                        dest="pca",
                        action='store_true',
                        help="Flag toggle use of PCA method ")

    parser.add_argument("-tsne",
                        dest="tsne",
                        action='store_true',
                        help="Flag toggle use of t-SNE method ")

    options = parser.parse_args()

    if not options.database:
        print(">>> A database location must be given.")
        parser.print_help()
        sys.exit(1)

    return options


def main():
    args = get_args()
    activities = get_chembl_activies(args.database, args.target_ids)

    target_freq = Counter(act.target_name for act in activities)
    for key, value in target_freq.items():
        print("{} activities found for {}".format(value, key))

    methods = []
    if args.tsne: methods.append("t-SNE")
    if args.pca: methods.append("PCA")

    cluster_reference_data(args.target_ids, activities, methods, args.out_file, args.input_smiles)


if __name__ == "__main__":
    main()
    print('\nMain() is done running.')
