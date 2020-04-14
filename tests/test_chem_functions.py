from unittest import TestCase
from MCR import chem_functions
from rdkit import Chem


class TestChemFunctions(TestCase):

    def test_remove_salts_1(self):
        example = ["CN(C)C.Cl", "CN(C)C"]

        result = chem_functions.remove_salts_mol(Chem.MolFromSmiles(example[0]))

        self.assertEqual(Chem.MolToSmiles(result), example[1])

    def test_remove_salts_2(self):
        input = "CC(=O)[O-].[Na+]"

        result = chem_functions.remove_salts_mol(Chem.MolFromSmiles(input))

        self.assertFalse(result.GetNumAtoms())

    def test_decharge_mol(self):
        test_mols = [Chem.MolFromSmiles(x) for x in ("CC(=O)[O-]", "c1ccc[n-]1")]

        result = [Chem.MolToSmiles(chem_functions.decharge_mol(x)) for x in test_mols]

        self.assertEqual(['CC(=O)O', 'c1cc[nH]c1'], result)

    def test_split_mol(self):
        test_mol = Chem.MolFromSmiles('CC(=O)[O-].[NH3+]C')

        result = chem_functions.split_mol(test_mol)

        self.assertTrue(len(result) == 2)

    def test_get_largest_fragment(self):
        test_mol = Chem.MolFromSmiles('CC(=O)[O-].[NH3+]C')

        largest_frag = chem_functions.get_largest_fragment_mol(test_mol)

        self.assertTrue('CC(=O)[O-]' == Chem.MolToSmiles(largest_frag, canonical=True))

    def test_substitute_reactive_group_1(self):
        target_mol = Chem.MolFromSmiles("CN(C)CC(Br)")
        frag = Chem.MolFromSmiles("Br")

        substituted_mol = chem_functions.substitute_reactive_group(target_mol, frag, 3)

        self.assertEqual(len(substituted_mol), 1)
        self.assertEqual(Chem.MolToSmiles(substituted_mol[0], canonical=True), "CN(C)CC[*:3]")

    def test_substitute_reactive_group_2(self):
        target_mol = Chem.MolFromSmiles("C1(=CC=C(C=C1C)C(=O)[H])C(=O)[H]")
        frag = Chem.MolFromSmarts("[CX3H]=O")
        result = set(['Cc1cc([*:3])ccc1C=O', "Cc1cc(C=O)ccc1[*:3]"])

        substituted_mol = chem_functions.substitute_reactive_group(target_mol, frag, 3)

        self.assertEqual(len(substituted_mol), 2)
        self.assertEqual((set([Chem.MolToSmiles(x, canonical=True) for x in substituted_mol])), result)

    def test_weld_r_groups_1(self):
        mol_to_weld = Chem.MolFromSmiles(
            "CN(C)CC(Br)c1cc([*:2])c([*:1])cn1.[H]C([H])([H])[*:1].[H][*:2]")

        welded_mol = chem_functions.weld_r_groups(mol_to_weld)

        self.assertEqual("Cc1ccc(C(Br)CN(C)C)nc1", Chem.MolToSmiles(welded_mol, canonical=True))

    def test_weld_r_groups_2(self):
        mol_to_weld = Chem.MolFromSmiles(
            "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C.C1=COC=C1[*:1].C[*:2]")

        welded_mol = chem_functions.weld_r_groups(mol_to_weld)

        self.assertEqual("CCC1=C(C)NC(=O)N[C@@H]1c1ccoc1", Chem.MolToSmiles(welded_mol, canonical=True))

    def test_generate_fingerprints(self):
        test_input = ['CC1=C(CC=O)[C@@H](C=O)NC(=O)N1', 'CCCC(O)CC1=C(C)NC(=O)N[C@@H]1C=O',
                      'CCC(=O)CC1=C(C)NC(=O)N[C@@H]1C=O']
        mols = [Chem.MolFromSmiles(x) for x in test_input]
        features = chem_functions.generate_fingerprints_iter(mols)

        # print(features.shape)
        # TOPOLISH

    def test_make_smiles_list(self):
        test_mol = Chem.MolFromSmiles('CC1=C(CC=O)[C@@H](C=O)NC(=O)N1')

        test_input = [[1, test_mol], ["hallo", [1, test_mol]]]

        self.assertEqual([[1, 'CC1=C(CC=O)[C@@H](C=O)NC(=O)N1'], ['hallo', [1, 'CC1=C(CC=O)[C@@H](C=O)NC(=O)N1']]],
                         chem_functions.make_smiles_list(test_input))

    def test_standardize_mol(self):
        target_mol = Chem.MolFromSmiles("CN(C)CC(Br)")

        result = Chem.MolToSmiles(chem_functions.standardize_mol(target_mol))

        self.assertEqual(result, 'CN(C)CCBr')

    def test_join_fragments(self):
        scaffold = "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C"
        sequence = [['CC(CC(C)(C)C)[*:1]', '000000'], ['O=C[*:2]', '000000']]

        result = chem_functions.join_fragments(sequence, scaffold)

        self.assertEqual(Chem.MolToSmiles(result[0]), "CC(CC(C)(C)C)[*:1].CC1=C(C[*:2])[C@H]([*:1])NC(=O)N1.O=C[*:2]")

    def test_lipinski(self):
        test_smiles = [
            'CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC',
            'CC1(CCC(=C(C1)c2ccc(cc2)Cl)CN3CCN(CC3)c4ccc(c(c4)'
            'c5cc6cc[nH]c6nc5)C(=O)NS(=O)(=O)c7ccc(c(c7)[N+](=O)[O-])NCC8CCOCC8)C'
        ]
        expected_results = [
            True,
            False
        ]

        test_mols = [Chem.MolFromSmiles(smiles) for smiles in test_smiles]
        result = [chem_functions.filter_lipinski(mol) for mol in test_mols]

        self.assertEqual(result, expected_results)

    def test_find_pains(self):
        from MCR.pains import pains_from_smarts
        pains = pains_from_smarts()
        test_smiles = [
            'C1=CC(=C(C=C1CCN)O)O',
            'O=C2C=CC(C1=CC=CC=C12)=O',
            'ClC=1C(=O)C(\Cl)=C(\Cl)C(=O)C=1Cl',
            'n2c1c(ncnc1n(c2)[C@@H]3O[C@@H]([C@@H](O)[C@H]3O)CO)N'
        ]
        expected_results = [
            False,
            False,
            False,
            True
        ]

        test_mols = [Chem.MolFromSmiles(smiles) for smiles in test_smiles]
        result = [chem_functions.filter_pains(test_mol, pains) for test_mol in test_mols]

        self.assertEqual(result, expected_results)
