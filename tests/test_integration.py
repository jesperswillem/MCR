import argparse
from unittest import TestCase
from unittest import mock  # python 3.3+

class TestIntegration(TestCase):
    def tearDown(self) -> None:
        import glob, shutil
        to_delete = glob.glob("../tests/delete_testdir*")
        for directory in to_delete:
            shutil.rmtree(directory)

    def test_washing_with_dask(self):
        """Bit more elaborate of a test to see if rdkit handles a set of molecules in a consistent way in coming versions
        and to see how dask handles the used chem_functions functions.
        """

        expected = ['CC(C)=CCC/C(C)=C\\CC/C(C)=C\\CO', 'CC12CC(O)C(CC1=O)C2(C)C', 'Oc1cc(C2CCNCC2)on1',
                    'Cn1ncc2cc(CN)ccc21', 'O=C(O)c1cc(Cl)cs1', 'Cc1cc(CN)ncc1Br', 'CO[C@@H](C)[C@@H](N)C(=O)O',
                    'Nc1ccc(Br)c(F)c1[N+](=O)[O-]', 'Cc1ccc(F)c(C#N)n1', 'Cc1ccc(F)c(CN)n1']

        from rdkit import Chem
        from MCR import chem_functions
        import dask.bag as db
        from rdkit import RDLogger, rdBase

        rdBase.DisableLog('rdApp.error')
        RDLogger.DisableLog('rdApp.info')

        bag = db.read_text("../tests/test_data/test_db.smi", blocksize=16e6)
        bag = bag.map(lambda x: Chem.MolFromSmiles(x)).filter(lambda x: x is not None)
        bag = bag.map(chem_functions.remove_salts_mol)
        bag = bag.map(chem_functions.decharge_mol)
        bag = bag.map(chem_functions.get_largest_fragment_mol)
        bag = bag.map(chem_functions.standardize_mol)

        self.assertEqual([Chem.MolToSmiles(x) for x in bag.take(10)], expected)

    def test_join_weld(self):
        from MCR import chem_functions
        from rdkit import Chem

        scaffold = "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C"
        sequence = [['CC(CC(C)(C)C)[*:1]', '000000'], ['O=C[*:2]', '000000']]

        result_mol = chem_functions.join_fragments(sequence, scaffold)
        result_smiles = Chem.MolToSmiles(chem_functions.weld_r_groups(result_mol[0]))

        with self.subTest(msg="Join and weld"):
            self.assertEqual(result_smiles, 'CC1=C(CC=O)[C@@H](C(C)CC(C)(C)C)NC(=O)N1')

