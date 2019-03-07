import sqlite3
from collections import namedtuple
from rdkit import Chem

from rdkit.Chem import SmilesMolSupplier, Fragments, SDWriter, AllChem, rdmolfiles, rdMolAlign, rdFMCS, rdchem, rdDistGeom, rdShapeHelpers


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


def create_images(mols, folder):
    for i, mol in enumerate(mols):
        Chem.Draw.MolToFile(mol, folder + "/norA_" + str(i) + ".png")


activity_tuple = namedtuple('activity_tuple', 'target_id, target_name, compound_id, canonical_smiles, standard_type, standard_relation, pchembl_value, standard_units')
def namedtuple_factory(cursor, row):
    return activity_tuple(*row)


def get_chembl_activies(db, targets):
    conn = sqlite3.connect(db)
    #setup factory so the rows are returned as named tuples.
    conn.row_factory = namedtuple_factory
    #load query for collecting bioactivity data.
    with open('/home/jespers/software/MCR/chembl_query.sql', 'r') as f:
        query = f.read()
    #execute query for each target.

    if type(targets) == str:
        return list(conn.execute(query, [targets]))

    result = []
    for t in targets:
        result += conn.execute(query, [t])
    return result


def is_int_or_float(a):
    try:
        float(a)
        return True
    except:
        return False


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
