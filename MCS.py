import os, argparse, sys
import dask.bag as db

from collections import defaultdict
from dask.diagnostics import ProgressBar
from rdkit import Chem, SimDivFilters, DataStructs, rdBase
from rdkit.Chem.rdchem import EditableMol

from shared import write_smi_file, create_images, get_smi_files



def filter_smiles(text, query_substruct_mol, max_h_atoms):
    mbag = text.map(lambda x: Chem.MolFromSmiles(x))
    mbag = mbag.filter(lambda x: x is not None).filter(lambda x: x.HasSubstructMatch(query_substruct_mol))
    mbag = mbag.filter(lambda x: x.GetNumHeavyAtoms() <= max_h_atoms)
    with ProgressBar():
        return mbag.compute()


def weld_r_groups(input_mol):
    # First pass loop over atoms and find the atoms with an AtomMapNum
    join_dict = defaultdict(list)
    for atom in input_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            join_dict[map_num].append(atom)

    # Second pass, transfer the atom maps to the neighbor atoms
    for idx, atom_list in join_dict.items():
        if len(atom_list) == 2:
            atm_1, atm_2 = atom_list
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()][0]
            nbr_1.SetAtomMapNum(idx)
            nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()][0]
            nbr_2.SetAtomMapNum(idx)

    # Nuke all of the dummy atoms
    new_mol = Chem.DeleteSubstructs(input_mol, Chem.MolFromSmarts('[#0]'))

    # Third pass - arrange the atoms with AtomMapNum, these will be connected
    bond_join_dict = defaultdict(list)
    for atom in new_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            bond_join_dict[map_num].append(atom.GetIdx())

    # Make an editable molecule and add bonds between atoms with correspoing AtomMapNum
    em = EditableMol(new_mol)
    for idx, atom_list in bond_join_dict.items():
        if len(atom_list) == 2:
            start_atm, end_atm = atom_list
            em.AddBond(start_atm, end_atm,
                       order=Chem.rdchem.BondType.SINGLE)

    final_mol = em.GetMol()

    # remove the AtomMapNum values
    for atom in final_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    final_mol = Chem.RemoveHs(final_mol)
    return final_mol


def substitute(compound, r_subs):
    mod_mol = []
    r=[[]]
    for x in r_subs:
        t = []
        for y in x:
            for i in r:
                t.append(i+[y])
        r = t
    sub_sets = []
    for i in t:
        sub_sets.append(".".join(i))
    for sub in sub_sets:
        mol_to_weld = Chem.MolFromSmiles(compound + '.' + sub)
        mod_mol.append(weld_r_groups(mol_to_weld))
    return mod_mol


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

    parser.add_argument("-e", "--enum_starting_point",
                        dest="enum_starting_point",
                        default=1,
                        nargs='*',
                        help="Queries for structures can be used as reactants, make sure these are in the right order.")

    options = parser.parse_args()

    if not options.database:
        print(">>> A database location must be given.")
        parser.print_help()
        sys.exit(1)
    elif not options.scaffold:
        print(">>> Scaffold structure is a required argument.")
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
    scaffold = args.scaffold

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
    write_smi_file(mols, 'outfile_result_mols.smi', args.enum_starting_point)

    if args.image_dir:
        print("writing images to {}".format(args.image_dir))
        if not os.path.isdir(args.image_dir):
            os.mkdir(args.image_dir)
        create_images(mols, args.image_dir)


if __name__ == "__main__":
    main()
    print('\nMain() is done running.')
