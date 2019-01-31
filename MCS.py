import os, shutil, gzip, sys, glob, argparse
import dask.bag as db

from dask.diagnostics import ProgressBar
from gen_ligands import *
from cluster_ligands import *

from io import StringIO
from operator import itemgetter

from rdkit import Chem, SimDivFilters, DataStructs, rdBase
from rdkit.Chem import SmilesMolSupplier, Fragments, SDWriter, AllChem, rdmolfiles, rdMolAlign, rdFMCS, rdchem, rdDistGeom, rdShapeHelpers
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

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
    with ProgressBar():
        return mbag.compute()

def get_smi_files(directory):
    if directory[-4:] == '.smi':
        return [directory]
    elif os.path.isdir(directory):
        files = []
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith('.smi'):
                    files.append(file)
        if len(files):
            return files
    print("A database location needs to either be a .smi file or directory containing .smi files")
    parser.print_help()
    sys.exit(1)

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("scaffold",
                        help="Smiles of the scaffold to perform MCS on.")

    parser.add_argument("-i", "--image_dir",
                        dest="image_dir",
                        default=None,
                        help="Directory to write images of found compounds to, will skip if left on default.")

    parser.add_argument("-o", "--out_file",
                        dest="out_file",
                        default="outfile.smi",
                        help="outfile to write found smiles to. Default is ./outfile.smi")

    parser.add_argument("-m", "--max_h_atoms",
                        dest="max_h_atoms",
                        default=50,
                        type=int,
                        help="Put an upper limit on number of heavy atoms in hits, default is 50.")

    parser.add_argument("-d", "--database",
                        dest="database",
                        default=None,
                        help="Directory to write images of found compounds to, will skip if left on default.")

    parser.add_argument("-q", "--queries",
                        dest="queries",
                        default=None,
                        nargs='*',
                        help="Queries for structures can be used as reactants, make sure these are in the right order.")

    options = parser.parse_args()

    if not options.database:
        print("A database location must be given.")
        parser.print_help()
        sys.exit(1)
    elif not options.scaffold:
        print("Scaffold structure is a required argument.")
        parser.print_help()
        sys.exit(1)

    return options


def main():
    # Read in arguments

    rdBase.DisableLog('rdApp.error')
    args = get_args()

    print("substituting wildcards in {} with:".format(args.scaffold))
    for i, q in enumerate(args.queries):
        print("{}. {}".format(i+1, q))
    print()

    # Read in database and apply specified filters using Dask for parallization/memory managment
    smiles_database = get_smi_files(args.database)
    text = db.read_text(smiles_database, blocksize =int(1e7))
    r_subs = []
    max_h_atoms = args.max_h_atoms

    for i, q in enumerate(args.queries):
        i += 1
        query_substruct_mol = Chem.MolFromSmarts(q)
        print("searching {} for compounds matching {}:".format(args.database, q))
        mols = filter_smiles(text, query_substruct_mol, max_h_atoms)
        print("found {}\n".format(len(mols)))
        write_smi_file(mols, 'outfile_fragment{}.smi'.format(i))
        mod_mols = []
        for mol in mols:
            # Generalize this for future use!!!
            mod_mol = Chem.ReplaceSubstructs(mol,
                                             query_substruct_mol,
                                             Chem.MolFromSmiles("[*:{}]".format(i)),
                                             )
            mod_mols.append(mod_mol[0])
        r_subs.append([Chem.MolToSmiles(x) for x in mod_mols])

    # Weld molecules
    print("Welding new mols...")
    mols = substitute(scaffold, r_subs)
    print("{} new mols generated, writing smiles to {}".format(len(mols), args.out_file))

    # Write result molecules to SD file
    write_smi_file(mols, 'outfile_result_mols.smi', enum_starting_point)

    if args.image_dir:
        print("writing images to {}".format(args.image_dir))
        if not os.path.isdir(args.image_dir):
            os.mkdir(args.image_dir)
        create_images(mols, args.image_dir)


if __name__ == "__main__":
    main()
    print("main() finished running without problems")
