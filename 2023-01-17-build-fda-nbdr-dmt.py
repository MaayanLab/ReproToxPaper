import pandas as pd
import requests
import itertools

### FDA ###
# extract FDA Category D and X drug data
fda_d = pd.read_csv('input/FDApreg_Category_D.tsv', sep='\t')
fda_x = pd.read_csv('input/FDApreg_Category_X.tsv', sep='\t')

# build DMT for FDA Cat D/X drugs
with open('data/FDA_Drugs.dmt', 'w') as f_out: 
  d = '\t'.join([term.strip().lower() for term in fda_d['Category D Drugs'].tolist()])
  x = '\t'.join([term.strip().lower() for term in fda_x['Category X Drugs'].tolist()])
  f_out.write(f"FDA Category D\t\t{d}\n")
  f_out.write(f"FDA Category X\t\t{x}")


### NBDR ###
# extract National Birth Defects Registry data
bd_otc = pd.read_csv('input/SBD_maternalExpo_BDorg_OTC.tsv', sep='\t', index_col=0).astype(int)
bd_pres = pd.read_csv('input/SBD_maternalExpo_BDorg_prescription.tsv', sep='\t', index_col=0, header=1).astype(int)

# fix formatting of drugs
nbd_drugs = bd_otc.index.tolist() + bd_pres.index.tolist()
nbd_drugs = [x.lower().strip().replace(' cream','') for x in nbd_drugs]
for i in range(len(nbd_drugs)): 
  if ' or ' in nbd_drugs[i]: 
    nbd_drugs.append(nbd_drugs[i].split('or')[1].strip())
    nbd_drugs[i] = nbd_drugs[i].split('or')[0].strip()

# query FDA for generic drug names
fda_endpoint = 'https://api.fda.gov/drug/label.json'
generic = {x: '' for x in nbd_drugs}

for x in generic.keys():
  resp = requests.get(f"{fda_endpoint}?search=openfda.brand_name:{x.lower()}").json()
  if 'error' in resp.keys():
    continue
  else:
    try:
      if len(resp['results']) > 0: 
        results = resp['results'][0]['openfda']['generic_name'][0].lower()
        generic[x] = [x.strip() for x in results.split(',')]
    except:
      print("No results for", x)

# get all NBDR drugs
generic_drugs = [x if generic[x] == '' else generic[x] for x in generic.keys()]
generic_drugs = set(itertools.chain.from_iterable(generic_drugs))

# built DMT for NBDR drugs
with open('data/NBDR_All_Drugs_Generic.dmt', 'w') as f_out: 
  f_out.write('NBDR\t\t' + '\t'.join(generic_drugs))