import os, shutil, gzip, sys, glob, argparse
import dask.bag as db
from dask.diagnostics import ProgressBar

from io import StringIO
from operator import itemgetter

from rdkit import Chem, SimDivFilters, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import rdBase

from shared import write_sd_file, write_smi_file, create_images


def filter_smiles_sim(text, query_substruct_mol, max_h_atoms, cutoff):
    query_fp = FingerprintMols.FingerprintMol(query_substruct_mol)
    mbag = text.map(lambda x: Chem.MolFromSmiles(x))
    mbag = mbag.filter(lambda x: x is not None)
    mbag = mbag.filter(lambda x: x.GetNumHeavyAtoms() <= max_h_atoms)
    mbag = mbag.filter(lambda x: DataStructs.FingerprintSimilarity(query_fp, FingerprintMols.FingerprintMol(x)) >= cutoff)
    with ProgressBar():
        return mbag.compute()


def filter_smiles_match(text, query_substruct_mol, max_h_atoms):
    mbag = text.map(lambda x: Chem.MolFromSmiles(x))
    mbag = mbag.filter(lambda x: x is not None).filter(lambda x: x.HasSubstructMatch(query_substruct_mol))
    mbag = mbag.filter(lambda x: x.GetNumHeavyAtoms() <= max_h_atoms)
    with ProgressBar():
        return mbag.compute()


def get_smi_files(d):
    if d[-4:] == '.smi':
        return [d]
    elif os.path.isdir(d):
        files = []
        for root, dirs, files in os.walk(d):
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
    
    parser.add_argument("query",
                        help  = "Smarts to perform query with.")

    parser.add_argument("-s", "--smilarity",
                        dest="similarity",
                        type=float,
                        default=None,
                        help="Supply a cut-off for similarity search based on 2d fingerprints, if not given will just match the substructures.")

    parser.add_argument("-i", "--image_dir",
                      dest    = "image_dir",
                      default = None,
                      help    = "Directory to write images of found compounds to, will skip if left on default.")
    
    parser.add_argument("-o", "--out_file",
                      dest    = "out_file",
                      default = "outfile.smi",
                      help    = "outfile to write found smiles to. Default is ./outfile.smi")
    
    parser.add_argument("-m", "--max_h_atoms",
                      dest    = "max_h_atoms",
                      default = 50,
                      help    = "Put an upper limit on number of heavy atoms in hits, default is 50.")
    
    parser.add_argument("-d", "--database",
                      dest    = "database",
                      default = None,
                      help    = "Directory to write images of found compounds to, will skip if left on default.")
    
    options = parser.parse_args()
    
    if not options.database:
        print(">>> A database location must be given.")
        parser.print_help()
        sys.exit(1)
    elif not options.query:
        print(">>> No smarts structure to query given.")
        parser.print_help()
        sys.exit(1)
    
    return options


def main():
    rdBase.DisableLog('rdApp.error')

    args = get_args()

    databases = [['/home/vonderent/databases/emolecules/emol_full.smi', False, "emol"]]
    query_sub = Chem.MolFromSmarts(args.query)

    smi_files = get_smi_files(args.database)
    text = db.read_text(smi_files, blocksize =int(1e7))

    print("searching {} for compounds matching {}:".format(args.database, args.query))

    if args.similarity:
        mols = filter_smiles_sim(text, query_sub, 80, args.similarity)
    else:
        mols = filter_smiles_match(text, query_sub, 80)

    print("found {} matches, writing results to drive...".format(len(mols)))

    if args.image_dir:
        if not os.path.isdir(args.image_dir):
            os.mkdir(args.image_dir)
        create_images(mols, args.image_dir)
    write_smi_file(mols, args.out_file)


if __name__ == "__main__":
    main()
    print('\nMain() is done running.')
