import argparse
import glob
import sys
from unittest import TestCase
from unittest import mock  # python 3.3+

sys.path.append('../')
import bin.mcr_run
import bin.mcr_reconstruct_clusters
import bin.mcr_plot


class TestFunction(TestCase):
    def setUp(self) -> None:
        self.to_delete_folders = []
        self.to_delete_files = []

    def tearDown(self) -> None:
        import shutil
        import os
        self.to_delete_folders += glob.glob("../tests/delete_testdir*")

        for directory in self.to_delete_folders:
            shutil.rmtree(directory)

        for f in self.to_delete_files:
            os.remove(f)

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="test_inputs/test1.inp"))
    def test_main1(self, mock_args):
        bin.mcr_run.main()

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="test_inputs/test2.inp"))
    def test_main2(self, mock_args):
        bin.mcr_run.main()

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="test_inputs/test3.inp"))
    def test_main3(self, mock_args):
        bin.mcr_run.main()

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(input_file="test_data/test_clusters/cluster.inp", clusters=[1, 2]))
    def test_main4(self, mock_args):
        bin.mcr_reconstruct_clusters.main()
        # self.to_delete_files += glob.glob('test_data/test_clusters/constructed_MCR_cluster_*.smi')
