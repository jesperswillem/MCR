import os, shutil, gzip, sys
import dask.bag as db

from gen_ligands import *
from cluster_ligands import *

from io import StringIO
from operator import itemgetter

from rdkit import Chem, SimDivFilters, DataStructs
from rdkit.Chem import SmilesMolSupplier, Fragments, SDWriter, AllChem, rdmolfiles, rdMolAlign, rdFMCS, rdchem, rdDistGeom, rdShapeHelpers
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from settings import *

def write_sd_file(mols, out_file):
    with open(out_file, 'w') as file:
        writer = SDWriter(file)
        for mol in mols:
            if mol == None: continue
            writer.write(Chem.MolToSmiles(mol))
        writer.close()


def write_smi_file(mols, out_file, start = 1):
    with open(out_file, 'w') as file:
        mols = filter(None, mols)
        for i, mol in enumerate(mols, start):
            file.write(Chem.MolToSmiles(mol) + ' ' + str(i).zfill(4) + '\n')

def filter_smiles(text, query_substruct_mol, max_h_atoms):
    mbag = text.map(lambda x: Chem.MolFromSmiles(x))
    mbag = mbag.filter(lambda x: x is not None).filter(lambda x: x.HasSubstructMatch(query_substruct_mol))
    mbag = mbag.filter(lambda x: x.GetNumHeavyAtoms() <= max_h_atoms)
    return mbag.compute()


def main():
    # Read in database and apply specified filters using Dask for parallization/memory managment
    r_subs = []
    text = db.read_text(smiles_database, blocksize =int(1e7))
    for i, q in enumerate(query):
        smarts, max_h_atoms = q
        i += 1
        query_substruct_mol = Chem.MolFromSmarts(smarts)
        mols = filter_smiles(text, query_substruct_mol, max_h_atoms)
        write_smi_file(mols, 'outfile_fragment{}.smi'.format(i))
        mod_mols = []
        for mol in mols:
            # Generelize this for future use!!!
            mod_mol = Chem.ReplaceSubstructs(mol,
                                             query_substruct_mol,
                                             Chem.MolFromSmiles("[*:{}]".format(i)),
                                             )
            mod_mols.append(mod_mol[0])
        r_subs.append([Chem.MolToSmiles(x) for x in mod_mols])
    # Weld molecules
    mols = substitute(scaffold, r_subs)
    # Write result molecules to SD file
    write_smi_file(mols, 'outfile_result_mols.smi', enum_starting_point)
    # generate Tsne plot
#    mol_plot_tsne(mols, 'tsne_plot.png')
#    mol_plot_pca(mols, 'pca_plot.png')


if __name__ == "__main__":
    main()
