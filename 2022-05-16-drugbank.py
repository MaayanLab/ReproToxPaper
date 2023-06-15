import xmltodict

with open('data/drugbank-full-database.xml', 'rb') as fr:
  D = xmltodict.parse(fr)

ensure_list = lambda L: L if type(L) == list else [L]
def only_id(i):
  if type(i) == str: return i
  elif type(i) == dict: return i['#text']
  elif type(i) == list: return only_id(i[0])

drug_enzyme_dmt = {}
drug_transporter_dmt = {}
drug_ids = {}

for d in D['drugbank']['drug']:
  drug_id = f"{d['name']} ({only_id(d['drugbank-id'])})"
  drug_ids[drug_id] = {
    d['name'],
    only_id(d['drugbank-id']),
    *(
      [
        f"{i['resource']}:{i['identifier']}"
        for i in ensure_list(d['external-identifiers']['external-identifier'])
      ] if d['external-identifiers'] else []
    ),
  }
  if d['enzymes']:
    drug_enzyme_dmt[drug_id] = {
      f"{e['name']} ({only_id(e['id'])})"
      for e in ensure_list(d['enzymes']['enzyme'])
    }
  if d['transporters']:
    drug_transporter_dmt[drug_id] = {
      f"{t['name']} ({only_id(t['id'])})"
      for t in ensure_list(d['transporters']['transporter'])
    }

def transpose_dmt(dmt):
  dm = {}
  for t, ds in dmt.items():
    for d in ds:
      if d not in dm:
        dm[d] = set()
      dm[d].add(t)
  return dm

enzyme_drug_dmt = {
  k: v
  for k, v in transpose_dmt(drug_enzyme_dmt).items()
  if len(v) >= 4
}
transporter_drug_dmt = {
  k: v
  for k, v in transpose_dmt(drug_transporter_dmt).items()
  if len(v) >= 4
}

with open('drugbank_enzyme.dmt', 'w') as fw:
  for term, ds in enzyme_drug_dmt.items():
    print(term, '', *ds, sep='\t', file=fw)

with open('drugbank_transporter.dmt', 'w') as fw:
  for term, ds in transporter_drug_dmt.items():
    print(term, '', *ds, sep='\t', file=fw)

with open('drugbank_drug_synonyms.dmt', 'w') as fw:
  for term, ds in drug_ids.items():
    print(term, '', *ds, sep='\t', file=fw)
