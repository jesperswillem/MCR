import os

# input database of for instance commercially available compounds
smiles_database = '/home/vonderent/databases/zinc/2018-11-27_ZINC_in_stock'
folder = True
#smiles_database = '/home/vonderent/2.lig_gen/1.A2b_ago/1.databases/emolecules/emol_full.smi'

# scaffold with numbered wildcards in postions to be 
scaffold = "C([*:2])C1[C@H]([*:1])NC(=O)NC=1C"
query = [["[CX3H1](=O)", 9],
         ["[CH2][CX3]([CH3])(=O)", 9]]

enum_starting_point = 13
