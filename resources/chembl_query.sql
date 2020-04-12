SELECT t.chembl_id AS target_chembl_id,
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
t.chembl_id = ?;