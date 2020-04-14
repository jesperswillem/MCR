import unittest, argparse
from unittest import mock, TestCase  # python 3.3+
import bin.mcr_run
from MCR import core


class TestCoreFunctions(unittest.TestCase):

    def test_parse_input_field_1(self):
        test_input = "[CX3H]=O, [CX3H](=O)N, C([*:2])C1[C@H]([*:1])NC(=O)NC=1C !Comment"
        desired_output = ["[CX3H]=O", "[CX3H](=O)N", "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C"]

        output = core.parse_input_field(test_input)

        self.assertEqual(desired_output, output)

    def test_parse_input_field_2(self):
        test_input = "[Br,I], N(O)=O"
        desired_output = ["[Br,I]", "N(O)=O"]

        output = core.parse_input_field(test_input)

        self.assertEqual(desired_output, output)

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="test_inputs/test1_non_rxn.inp"))
    def test_parse_input(self, mock_args):
        settings, reactants, MCR_parameters = bin.mcr_run.parse_input()

        self.assertEqual(settings["num_processes"], 4)
        self.assertEqual(reactants[0]["include_smarts"], "")
        self.assertEqual(len(reactants[0]["exclude_smarts"]), 2)
        self.assertFalse(MCR_parameters['perform_mcr'])

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="does/not/exist/for/sure"))
    def test_parse_input_sanityCheck(self, mock_args):
        with self.assertRaises(FileNotFoundError):
            bin.mcr_run.parse_input()

    def test_molbag_to_featurearray(self):
        import dask.bag as db

        with self.subTest(status_code='empty molbag') and self.assertRaises(ValueError):
            bag = db.from_sequence([])
            core.molbag_to_featurearray(bag)

        with self.subTest(status_code='flawed molbag') and self.assertRaises(Exception):
            bag = db.from_sequence([False, 1, 'mol'])
            core.molbag_to_featurearray(bag)

        with self.subTest(status_code='proper molecules'):
            from rdkit import Chem
            bag = db.read_text("../tests/test_data/aldehyde_query.smi", blocksize=16e6)
            mol_bag = bag.map(lambda x: Chem.MolFromSmiles(x)).filter(lambda x: bool(x))

            feature_array = core.molbag_to_featurearray(mol_bag)
            self.assertTrue((88, 512), feature_array.shape)

    def test_cull_empty_partitions(self):
        import dask.bag as db
        test_input = db.from_sequence([None, 1, 'test', [], 'test'])
        test_input = test_input.filter(lambda x: bool(x))

        result = core.cull_empty_partitions(test_input)

        self.assertEqual(3, result.npartitions)

    def test_construct_query(self):
        import dask.bag as db
        from rdkit import Chem

        input_file = 'test_data/aldehyde_query.smi'
        input_bag = db.read_text(input_file).map(lambda x: Chem.MolFromSmiles(x)).filter(lambda x: x != None)

        result = []
        for i in [5, 8, 12]:
            query = {
                "reacting_group": None,
                "include_smarts": None,
                "exclude_smarts": None,
                "max_heavy_ratoms": i,
                "wash_molecules": True,
                "keep_largest_fragment": True
            }
            bag = core.construct_query(query, input_bag)
            result.append(bag.count().compute())

        self.assertEqual(result, [0, 17, 57])

    def test_exclude_smarts(self):
        import dask.bag as db
        from rdkit import Chem

        input_file = 'test_data/aldehyde_query.smi'
        input_bag = db.read_text(input_file).map(lambda x: Chem.MolFromSmiles(x)).filter(lambda x: x != None)

        result = []

        query = {
            "reacting_group": None,
            "include_smarts": None,
            "exclude_smarts": ['[CX3H](O)[N;n]', '[Br,I]', '[CX2]#[NX1]'],
            "max_heavy_ratoms": 10,
            "wash_molecules": True,
            "keep_largest_fragment": True
        }
        bag = core.construct_query(query, input_bag)
        bag = bag.map(lambda x: Chem.MolToSmiles(x))

        for smiles in bag.compute():
            self.assertNotIn('Br', smiles)

    def test_exclude_smarts2(self):
        import dask.bag as db
        from rdkit import Chem

        input_file = 'test_data/aldehyde_query.smi'
        input_bag = db.read_text(input_file).map(lambda x: Chem.MolFromSmiles(x)).filter(lambda x: x != None)

        result = []

        query = {
            "reacting_group": None,
            "include_smarts": None,
            "exclude_smarts": '[Br,I]',
            "max_heavy_ratoms": 10,
            "wash_molecules": True,
            "keep_largest_fragment": True
        }
        bag = core.construct_query(query, input_bag)
        bag = bag.map(lambda x: Chem.MolToSmiles(x))

        for smiles in bag.compute():
            self.assertNotIn('Br', smiles)


    def test_execute_mcr(self):
        from rdkit import Chem
        reactants_list = [
            [('O=CC', 1), ('O=CBr', 2)],
            [('CC(C=O)CC(C)(C)C', 1)]
        ]
        reactants_list = [
            [(Chem.MolFromSmiles(a), b) for a, b in x] for x in reactants_list
        ]

        scaffold = ['C([*:2])C1[C@H]([*:1])NC(=O)NC=1C']
        reacting_groups = ['[CX3H]=O', '[CH2][CX3]([CH3])(=O)']

        print(core.execute_mcr(reactants_list, reacting_groups, scaffolds=scaffold)[0].compute())

    def test_execute_mcr_rxn(self):
        from rdkit import Chem
        from MCR import chem_functions, helpers
        test_reactants = [
            [('CC(C)OC(=O)C(C)(=O)', 1)],
            [('c1cocc1C(=O)', 1), ('c1ccsc1C(=O)', 2), ('O=CCC(F)C=O', 3)],
            [('NC(=O)N', 1), ('NC(=S)N', 2)]
        ]
        expected_results = [
            ['CC1=C(OC(C)C)[C@H](c2ccoc2)NC(=O)N1', [1, 1, 1]],
            ['CC1=C(OC(C)C)[C@H](c2ccoc2)NC(=S)N1', [1, 1, 2]],
            ['CC1=C(OC(C)C)[C@H](c2cccs2)NC(=O)N1', [1, 2, 1]],
            ['CC1=C(OC(C)C)[C@H](c2cccs2)NC(=S)N1', [1, 2, 2]],
            ['CC1=C(OC(C)C)[C@H](CC(F)C=O)NC(=O)N1', [1, 3, 1]],
            ['CC1=C(OC(C)C)[C@H](C(F)CC=O)NC(=O)N1', [1, 3, 1]],
            ['CC1=C(OC(C)C)[C@H](CC(F)C=O)NC(=S)N1', [1, 3, 2]],
            ['CC1=C(OC(C)C)[C@H](C(F)CC=O)NC(=S)N1', [1, 3, 2]]
        ]


        test_reactants = [
            [[Chem.MolFromSmiles(a), b] for a, b in x] for x in test_reactants
        ]

        scaffold = '[*:2]C1[C@H]([#6:1])[N]C(=[O,S:3])[N]C=1C'
        reactants = ['O=[CX3]([CH3])C[*:2]', '[#6:1][CX3H1]=O', '[N]C(=[O,S:3])[N]']
        # reactants = '.'.join(reactants)
        # reaction = f"{reactants}>>{scaffold}"

        unformatted_results = core.execute_mcr_rxn(test_reactants, reactants, [scaffold])[0].compute()
        results = [[Chem.MolToSmiles(i[0]), i[1]] for i in unformatted_results]

        self.assertEqual(expected_results, results)

