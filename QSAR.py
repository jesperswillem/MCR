import sys, csv, argparse
from collections import defaultdict, namedtuple

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Chem.rdchem import EditableMol

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import r2_score, mean_squared_error, matthews_corrcoef
from sklearn.model_selection import train_test_split, KFold, cross_val_score, cross_val_predict
from rdkit.ML.Scoring.Scoring import CalcBEDROC, CalcEnrichment

from shared import get_chembl_activies

def score_matt(model, x, y):
    pred = model.predict(x)
    matt = matthews_corrcoef(y, pred)
    return matt


def score_bedroc(model, x, y):
    prob = model.predict_proba(x)
    prob = [x[1] for x in prob]
    y_ord = [[y] for _, y in sorted(zip(prob, y), key=lambda pair: pair[0], reverse=True)]
    bedroc = CalcBEDROC(y_ord, 0, 20)
    return bedroc


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--out_file",
                        dest="out_file",
                        default=None,
                        help="outfile to write smiles and predictions to.")

    parser.add_argument("-i", "--input_smiles",
                        dest = "input_smiles",
                        default = None,
                        help = "Input smiles for which to make predictions.")

    parser.add_argument("-db", "--database",
                        dest="database",
                        default=None,
                        help="path to SQLite .db file of the ChEMBL database.")

    parser.add_argument("-t", "--target_id",
                        dest="target_id",
                        default=None,
                        help="ChEMBL id of the target for which you want to create a QSAR model.")

    options = parser.parse_args()

    if not options.database:
        print(">>> A database location must be given.")
        parser.print_help()
        sys.exit(1)
    elif options.input_smiles:
        if not options.out_file:
            print(">>> No outfile defined for compound predictions")
            parser.print_help()
            sys.exit(1)
    return options


def main():
    args = get_args()

    print("\nLoading chembl data...")
    chembl_data = get_chembl_activies(args.database, args.target_id)
    print("Found {} compounds for {}".format(len(chembl_data), args.target_id))

    print('\ncalculating fingerprints')
    mols = [Chem.MolFromSmiles(x.canonical_smiles) for x in chembl_data]
    X = [AllChem.GetMorganFingerprintAsBitVect(m, radius = 3, nBits = 256) for m in mols]

    Y_continuous = [x.pchembl_value for x in chembl_data]
    Y = [1 if x > 6.5 else 0 for x in Y_continuous]

    print("\nTraining QSAR random forest classifier...")
    rf = RandomForestClassifier(n_estimators=1000, random_state=8121993, n_jobs=16, verbose=0, max_depth=7)
    rf.fit(X, Y)

    print('\nPerforming 5-fold cross-validation...')
    scores = cross_val_score(rf, X, Y, cv=5, scoring=score_matt)
    print("Mattews: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    scores = cross_val_score(rf, X, Y, cv=5, scoring=score_bedroc)
    print("Bedroc: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

    if args.input_smiles:
        if not args.out_file:
            print("No outfile defined!")
            sys.exit(1)
        with open(args.input_smiles) as file:
            gen_mols = []
            for row in file:
                row = row[:-1]
                gen_mols.append(row.split(" "))

        mols = [Chem.MolFromSmiles(x[0]) for x in gen_mols]
        X_new = [AllChem.GetMorganFingerprintAsBitVect(m, radius=3, nBits=256) for m in mols]
        Y_pred_new = rf.predict_proba(X_new)
        gen_mols = [[x[0], x[1], Y_pred_new[i][0]] for i, x in enumerate(gen_mols)]
        with open(args.out_file, "w") as file:
            writer = csv.writer(file)
            writer.writerow(["canonical_smiles", "id", "predicted_proba"])
            for line in gen_mols:
                writer.writerow(line)

    """
    To be restored to functional later.
    print("predicting for test set")
    Y_pred = cross_val_predict(rf, X, Y, cv=5, method='predict_proba')
    print(Y_pred)
    Y_pred = [x[1] for x in list(Y_pred)]
    print(chembl_data[:5])
    chembl_output = [[x.canonical_smiles, x.compound_id, Y_pred[i], x[7], x[2]] for i, x in enumerate(chembl_data)]
    with open("chembl_mol_predictions.smi", "w") as file:
        writer = csv.writer(file)
        writer.writerow(["canonical_smiles", "id", "predicted_proba", "true_bin", "p_ChEMBL"])
        for line in chembl_output:
            writer.writerow(line)

    print("loading generated compounds")
    with open(args.input_smiles) as file:
        gen_mols = []
        for row in file:
            row = row[:-1]
            gen_mols.append(row.split(" "))
    mols = [Chem.MolFromSmiles(x[0]) for x in gen_mols]
    X_new = [AllChem.GetMorganFingerprintAsBitVect(m, radius=3, nBits=256) for m in mols]
    rf.fit(X, Y)
    Y_pred_new = rf.predict_proba(X_new)
    Y_pred_new = [x[0] for x in Y_pred_new]
    gen_mols = [[x[0], x[1], Y_pred_new[i]] for i, x in enumerate(gen_mols)]
    with open("gen_mol_predictions.smi", "w") as file:
        writer = csv.writer(file)
        writer.writerow(["canonical_smiles", "id", "predicted_proba"])
        for line in gen_mols:
            writer.writerow(line)
    """

if __name__ == "__main__":
    main()
    print('\nMain() is done running.')
