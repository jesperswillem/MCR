"""
MCR chem_functions

Authors: Florian van der Ent, Willem Jespers
"""
from collections import defaultdict

import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.SaltRemover import SaltRemover



# decharge reactions taken from rdkit manual
decharge_reactions = [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )]

remover = SaltRemover()


def make_smiles_list(l):
    if type(l) == Chem.rdchem.Mol:
        return Chem.MolToSmiles(l)
    elif type(l) == list:
        return [make_smiles_list(j) for j in l]
    else:
        return l


def standardize_mol(mol):
    """Accepts `mol` rdkit molecule returns a sanitized version with properly assigend stereochemistry.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol object
    """
    flag1 = Chem.SanitizeMol(mol)
    flag2 = Chem.AssignStereochemistry(mol)
    return mol


def decharge_mol(mol):
    """Accepts `mol` rdkit molecule returns a decharged version if one of the smart strings is recognized. Based on code
    from the rdkit manual. 

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol object
    """
    if decharge_reactions is None:
        raise ValueError

    for i, (reactant, product) in enumerate(decharge_reactions):
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol


def split_mol(mol):
    """Accepts `mol` rdkit molecule returns any unconnected fragments as separate mol objects.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object

    Returns
    -------
    tuple : tuple[rdkit mol]
    """
    return Chem.rdmolops.GetMolFrags(mol, asMols=True)


def remove_salts_mol(mol):
    """Accepts `mol` rdkit molecule returns fragment with highest num heavy atoms.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object

    Returns
    -------
    largest fragment : rdkit.Chem.rdchem.Mol object or None
    """
    result = remover.StripMol(mol)

    return result


def custom_remove_salts_mol(mol):
    # split mol into disconnected fragments
    frags = split_mol(mol)

    # Remove single atom fragments and some common bigger ions.
    result = []
    if len(frags) == 1: return mol
    for frag in frags:
        frag_canon_smiles = Chem.MolToSmiles(frag)
        if frag_canon_smiles in ["[K+]",
                                 "[Na+]",
                                 "[Cl]",
                                 "[Br]",
                                 "[Ca2+]"]:
            continue
        if frag_canon_smiles in ['O=N(O)O',
                                 'O=P(O)(O)O',
                                 'F[PH](F)(F)(F)(F)F',
                                 'O=S(=O)(O)O',
                                 'CS(=O)(=O)O',
                                 'Cc1ccc(S(=O)(=O)O)cc1',
                                 'CC(=O)[O-]',
                                 'O=C(O)C(F)(F)F',
                                 'O=C(O)C=CC(=O)O',
                                 'O=C(O)C(=O)O',
                                 'O=C(O)C(O)C(O)C(=O)O',
                                 'C1CCC(NC2CCCCC2)CC1']:
            continue
        result.append(frag)


def get_largest_fragment_mol(mol):
    """Accepts `mol` rdkit molecule returns fragment with highest num heavy atoms.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object

    Returns
    -------
    largest fragment : rdkit.Chem.rdchem.Mol object
    """
    frags = split_mol(mol)
    if len(frags) == 1: return mol
    return sorted(frags, key = lambda x: x.GetNumHeavyAtoms(), reverse=True)[0]


def join_fragments(sequence, scaffold):
    """ Takes a sequence of [[mol, int], [mol, int], ...] and creates a single of the mols contained combined with the
    scaffold.

    Parameters
    ----------
    sequence : list(list), : rdkit.Chem.rdchem.Mol objectd

    Returns
    -------
    result : list(mol, int, int, etc)
    """

    smiles = [scaffold]
    result = []

    for frag, key in sequence:
        result.append(key)
        smiles.append(frag)

    mol = Chem.MolFromSmiles(".".join(smiles))
    result.insert(0, mol)
    return result


def substitute_reactive_group(target_mol, reactive_group_mol, number):
    """Accepts 2 `mol` rdkit molecule objects one target molecule and a reactive group which is present in the target
    molecule and replaces the reactive group in the target molecule with a numbered wildcard.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object, mol : rdkit.Chem.rdchem.Mol object, number : int

    Returns
    -------
    result : tuple
    """
    result = Chem.ReplaceSubstructs(
        target_mol,
        reactive_group_mol,
        Chem.MolFromSmiles(f"[*:{number}]")
    )

    return result


def generate_fingerprints_iter(mol_iter, nBits=512):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 4, nBits=nBits) for m in mol_iter]
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    return np.stack(np_fps)


def generate_fingerprints_iter_debug(mol_iter, nBits=512):
    fps = []
    for mol in mol_iter:
        try:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=nBits))
        except:
            print(Chem.MolToSmiles(mol))
            Chem.GetSymmSSSR(mol)
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=nBits))
    try:
        np_fps = []
        for fp in fps:
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            np_fps.append(arr)
        return np.stack(np_fps)
    except:
        return []


def generate_fingerprints(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=512)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def apply_rxn(reactants, rxn):
    """takes an rxn reaction string and a list of lists of reactants(rdkit mol objects) and returns a list of products.

    Parameters
    ----------
    reactans: list(str)
    rxn: str

    Returns
    -------
    reaction: list(list)
    """
    reaction = AllChem.ReactionFromSmarts(rxn)
    results = []
    for reactant in reactants:
        unfiltered_products = reaction.RunReactants(reactant[0],)

        found_smiles = set()
        products = []
        for product in unfiltered_products:
            product = product[0]
            smiles = Chem.MolToSmiles(product)
            if not smiles in found_smiles:
                # Chem.GetSSSR(product)
                Chem.SanitizeMol(product)
                products.append(product)
                found_smiles.add(smiles)

        results.append([(product, reactant[1]) for product in products])
    return results


def create_reaction(scaffold, *reactants):
    """Creates an rdkit reaction from a product scaffold and list of reactant smarts.

    Parameters
    ----------
    scaffold: str
    reactans: list(str)

    Returns
    -------
    reaction: rdkit.reaction
    """
    reactants = '.'.join(reactants)
    reaction = f"{reactants}>>{scaffold}"
    return AllChem.ReactionFromSmarts(reaction,)


def weld_r_groups(mol):
    """Accepts `mol` rdkit molecule of a scaffold and with numbered wildcards and fragments with corresponding
    wildcards and returns a fused molecule. Based on code posted by Patrick Walters on the RDkit forums.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object

    Returns
    -------
    final_mol : rdkit.Chem.rdchem.Mol object
    """

    # First pass loop over atoms and find the atoms with an AtomMapNum
    join_dict = defaultdict(list)
    for atom in mol.GetAtoms():
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
    new_mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('[#0]'))

    # Third pass - arrange the atoms with AtomMapNum, these will be connected
    bond_join_dict = defaultdict(list)
    for atom in new_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            bond_join_dict[map_num].append(atom.GetIdx())

    # Make an editable molecule and add bonds between atoms with corresponding AtomMapNum
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

def filter_lipinski(mol):
    """Takes a rdkit mol object and returns true if it conforms to the lipinski rules.

    Parameters
    ----------
    mol: rdkit:mol

    Returns
    -------
    boolean
    """
    if Lipinski.NumHAcceptors(mol) > 10:
        return False
    elif Lipinski.NumHDonors(mol) > 5:
        return False
    elif Descriptors.MolWt(mol) > 500:
        return False
    elif Descriptors.MolLogP(mol) > 5:
        return False

    return True


def filter_pains(mol, pains):
    """Takes a rdkit mol object and returns true if it contains no pains.

    Parameters
    ----------
    mol: rdkit:mol

    Returns
    -------
    boolean
    """
    hits = []
    for name, pain in pains:
        if mol.HasSubstructMatch(pain):
            hits.append(name)
    if hits:
        # print(hits)
        return False
    else:
        return True
