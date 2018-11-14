import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA



def mol_plot_pca(mols, path):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 4) for m in mols]
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    X = PCA(n_components=2).fit_transform(np_fps)
    plt.scatter(X[:, 0], X[:, 1], c="r")
    plt.savefig(path)

def mol_plot_tsne(mols, path):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 4) for m in mols]
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    X = TSNE(n_components=2).fit_transform(np_fps)
    plt.scatter(X[:, 0], X[:, 1], c="r")
    plt.savefig(path)

def mol_plot_tsne_subset(mols, path, subset):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 4) for m in mols]
    np_subset = []
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    subset = [subset[x] for x in range(len(subset))]
    X = TSNE(n_components=2, verbose=1).fit_transform(np_fps)
    plt.scatter(X[:, 0], X[:, 1], c="r")
    plt.savefig(path)

if __name__ == "__main__":
    pass
