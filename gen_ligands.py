import sys

import dask.bag as db

from itertools import product
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from collections import defaultdict
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Chem.rdchem import EditableMol


def remove_duplicates(mols):
    seen = set()
    num_in = len(mols)
    result_mols = []
    for mol in mols:
        can_smiles = Chem.MolToSmiles(mol, canonical = True)
        if can_smiles in seen:
            continue
        else:
            seen.add(mol)
            result_mols.append(mol)
    print('removed %s duplicate mols from list.' % str(num_in - len(result_mols)))
    return result_mols

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


if __name__ == "__main__":
    scaffold = 'C([*:2])1C([*:1])NC(O)NC=1([*:3])'
    r_subs = [["[*:1]c2ccco2", "[*:1]c2cccs2", "[*:1]C2CCCC2"],
              ["CC(C)OC(O)[*:2]"],
              ["C(F)(F)(F)[*:3]"]]
    mod_mol = substitute(scaffold, r_subs)
    for i, mol in enumerate(mod_mol):
        print(mol)
        print(Chem.MolToSmiles(mol))
        Draw.MolToFile(mol,'images/19_sub' + str(i) + '.png')
