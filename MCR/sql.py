from collections import namedtuple
import sqlite3


chembl_query = """SELECT t.chembl_id AS target_chembl_id,
t.pref_name AS target_name,
m.chembl_id AS compound_chembl_id,
s.canonical_smiles,
c.standard_type,
c.standard_relation,
c.pchembl_value,
c.standard_units
FROM target_dictionary t,
assays a,
activities c,
compound_structures s,
molecule_dictionary m,
compound_properties p
WHERE t.tid = a.tid AND
a.assay_id = c.assay_id AND
c.molregno = m.molregno AND
m.molregno = s.molregno AND
m.molregno = p.molregno AND
a.confidence_score IN (7, 9) AND
t.target_type = 'SINGLE PROTEIN' AND
a.assay_organism = "Homo sapiens" AND
c.pchembl_value IS NOT NULL AND
c.standard_relation = '=' AND
p.full_mwt < 900 AND
t.chembl_id = ?;"""


def get_chembl_activies(db, targets):

    activity_tuple = namedtuple('activity_tuple',
                                'target_id, target_name, compound_id, canonical_smiles, standard_type, '
                                'standard_relation, pchembl_value, standard_units')

    def namedtuple_factory(cursor, row):
        return activity_tuple(*row)

    conn = sqlite3.connect(db)
    # setup factory so the rows are returned as named tuples.
    conn.row_factory = namedtuple_factory
    # load query for collecting bioactivity data.
    # with open(f'{realpath}../resources/chembl_query.sql') as f:
    #     query = f.read()
    # execute query for each target.

    if type(targets) == str:
        return list(conn.execute(chembl_query, [targets]))

    result = []
    for t in targets:
        result += conn.execute(chembl_query, [t])
    return result
