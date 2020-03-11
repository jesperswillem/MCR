import unittest, argparse
from unittest import mock, TestCase  # python 3.3+
import bin.run_mcr
from MCR import core


class TestCoreFunctions(unittest.TestCase):

    def test_parse_input_field_1(self):
        test_input = "[CX3H]=O, [CX3H](=O)N, C([*:2])C1[C@H]([*:1])NC(=O)NC=1C !Comment"
        desired_output = ["[CX3H]=O", "[CX3H](=O)N", "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C"]

        output = core.parse_input_field(test_input)

        self.assertEqual(desired_output, output)

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="test_inputs/test1.inp"))
    def test_parse_input(self, mock_args):
        settings, reactants, MCR_parameters = bin.run_mcr.parse_input()

        self.assertEqual(settings["num_processes"], 4)
        self.assertEqual(reactants[0]["include_smarts"], "")
        self.assertEqual(len(reactants[0]["exclude_smarts"]), 2)
        self.assertFalse(MCR_parameters['perform_mcr'])

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="does/not/exist/for/sure"))
    def test_parse_input_sanityCheck(self, mock_args):
        with self.assertRaises(FileNotFoundError):
            bin.run_mcr.parse_input()

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
            bag = db.read_text("/Users/florian/Documents/uppMCR/tests/test_data/aldehyde_query.smi", blocksize=16e6)
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

