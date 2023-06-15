#%%
import re
import json
import uuid
import pathlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
from matplotlib import pyplot as plt
from supervenn import supervenn

#%%
class IDMapper:
    ''' Stores id mappings and makes it easy to use many of them in tandem
    '''
    def __init__(self):
        # { uuid1: {id1, id2, ...} }
        self._forward = {}
        # { id1: uuid1, id2, uuid1, ... }
        self._reverse = {}
        # { uuid1: { ns1: id1 }, ... }
        self._namespaces = {}
        # { ns1: { shared_synonym: { conflictid1: origid1 }, ... } } }
        self._conflicts = {}

    def summary(self):
        return Counter(
            frozenset(ns_ids.keys())
            for ns_ids in self._namespaces.values()
        )
    
    def conflicts_summary(self):
        return Counter({
            k: len(v)
            for k, v in self._conflicts.items()
        })

    def top_conflicts(self, k=10):
        return Counter({
            f"{ns}: {conflict}": len(cases)
            for ns, cc in self._conflicts.items()
            for conflict, cases in cc.items()
        }).most_common()[:k]

    
    def get_id(self, id, namespace=None):
        if id is None: return None
        if namespace is None:
            return dict(
                id=id,
                refs=self._namespaces[id],
                synonyms=self._forward[id],
            )
        else:
            return self._namespaces[id].get(namespace)

    def get(self, term, namespace=None):
        id = self._reverse.get(term)
        return self.get_id(id, namespace=namespace)

    def find(self, term):
        potential_ids = {
            id
            for k, id in self._reverse.items()
            if str(term).lower().strip() in str(k).lower().strip() or str(k).lower().strip() in str(term).lower().strip()
        }
        return {
            id: self.get_id(id)
            for id in potential_ids
        }

    def update(self, mappings, namespace=None):
        ''' Add mappings of the form:
        { identifier: { synonyms } }
        '''
        for key, synonyms in (mappings.items() if type(mappings) == dict else mappings):
            id = uuid.uuid4()
            self._forward[id] = Counter()
            self._namespaces[id] = {namespace: key}
            for synonym in {key, *synonyms}:
                if synonym not in self._reverse:
                    self._forward[id].update([synonym])
                    self._reverse[synonym] = id
                else:
                    orig_id = self._reverse[synonym]
                    if orig_id == id:
                        self._forward[id].update([synonym])
                    else:
                        for ns, k in self._namespaces.pop(id, {}).items():
                            if orig_id not in self._namespaces:
                                self._namespaces[orig_id] = {}
                            orig_k = self._namespaces[orig_id].get(ns)
                            if orig_k is not None:
                                if orig_k != k:
                                    if ns not in self._conflicts:
                                        self._conflicts[ns] = {}
                                    if synonym not in self._conflicts[ns]:
                                        self._conflicts[ns][synonym] = {}
                                    self._conflicts[ns][synonym][k] = orig_k
                            else:
                                self._namespaces[orig_id][ns] = k
                        new_cnt = self._forward.pop(id)
                        self._forward[orig_id] += new_cnt
                        self._reverse.update({s: orig_id for s in new_cnt.keys()})
                        id = orig_id

def read_gmt_iter(p: pathlib.Path):
    with pathlib.Path(p).open('r') as fr:
        for line in fr:
            term, *geneset = line.strip('\n').split('\t')
            yield term, set(filter(None, geneset))

def read_gmt_dict(p: pathlib.Path):
    return dict(read_gmt_iter(p))

def one(it):
    try:
        return next(iter(it))
    except TypeError:
        return None
    except StopIteration:
        return None


def rdkit_canonical_smiles(smiles):
    from rdkit import Chem, RDLogger
    if not smiles or smiles == '-': return None
    try:
        RDLogger.DisableLog('rdApp.*')
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    except KeyboardInterrupt:
        raise
    except:
        return smiles

# %%
df_D_X = pd.read_excel('input/pregnancy_category_D_and_X_v2.xlsx')
df_D_X

#%%
df_drugs_com = pd.read_csv('data/drugs-com.tsv', sep='\t')
df_drugs_com['drug'] = df_drugs_com.drug.str.split('/')
df_drugs_com = df_drugs_com.explode('drug').reset_index(drop=True)
df_drugs_com['drug'] = df_drugs_com['drug'].str.strip()
df_drugs_com['drug'] = df_drugs_com['drug'].str.replace(' topical', '')
df_drugs_com

#%%
df_placenta = pd.read_excel('input/Placenta_Barrier.xlsx')
df_placenta

# %%
coexpression_matrix = np.load('data/L1000_2021_drug_similarity.npz', allow_pickle=True)
coexpression_matrix_index = coexpression_matrix['index'][:]

# %%
lincs_drug_metadata = pd.read_csv('data/Drugs_metadata.csv', index_col=0)
lincs_sm_meta = pd.read_csv('data/LINCS_small_molecules.tsv', sep='\t', index_col=0)
drugbank_synonyms = read_gmt_dict('data/drugbank_drug_synonyms.dmt')
with open('data/chembl_drugs.json', 'r') as fr:
    chembl_drugs = json.load(fr)

# %%
mapper = IDMapper()

def filter_none(it):
    for el in it:
        if (
            el is not None
            and (type(el) == str and el != '' and el != '(null)' and el != '-')
            and not pd.isnull(el)
        ):
            yield el

name_paren_extract = re.compile(r"^(.+?)( \(([^\)]+)\))?$")
def extract_names(s):
    if type(s) != str or not s: return []
    m = name_paren_extract.match(s)
    return filter_none([
        s,
        m and m.group(1),
        m and m.group(3),
    ])

mapper.update({
  drug['chemblId']: set(filter_none({
    f"chembl:{drug['chemblId'].lower()}",
    drug['prefName'].lower()
  } | {synonym.lower() for synonym in drug['synonyms']}))
  for drug in chembl_drugs
}, namespace='chembl')

mapper.update({
    k: set(filter_none(set.union(
        {v.strip().lower() for v in V},
        {k.strip().lower()},
        {n for v in V for n in extract_names(v.strip().lower())},
    ))) - {'drugs product database (dpd):331'}
    for k, V in drugbank_synonyms.items()
}, namespace='drugbank')

mapper.update({
    index.strip().lower(): set(filter_none({
        row['pert_name'].strip().lower(),
        rdkit_canonical_smiles(row['canonical_smiles']),
        row['inchi_key'],
    }))
    for index, row in lincs_sm_meta.iterrows()
}, namespace='lincs_sm_meta')

mapper.update({
    index.strip().lower(): set(filter_none({
        row['LSM_id'].strip().lower() if type(row['LSM_id']) == str else None,
        row['pert_iname'].strip().lower() if type(row['pert_iname']) == str else None,
        *extract_names(row['alt_name'].strip().lower() if type(row['alt_name']) == str else None),
        f"cid:{row['pubchem_cid']}" if type(row['pubchem_cid']) == str else None,
        rdkit_canonical_smiles(row['canonical_smiles']),
    }))
    for index, row in lincs_drug_metadata.iterrows()
}, namespace='lincs_drug_meta')

mapper.update({
    idx: set(filter_none({
        idx.strip().lower(),
        idx.replace('-',' ').strip().lower(),
    }))
    for idx in coexpression_matrix_index
}, namespace='lincs_matrix')

mapper.update({
    k: set(filter_none(v))
    for k, v in read_gmt_dict('data/2022-08-29-drug_synonmys.dmt').items()
}, namespace='structure')

mapper.update({
    row['name'].strip().lower(): set(filter_none({
        row['name'].strip().lower(),
        f"wikipedia:{row['name'].strip().lower()}",
        rdkit_canonical_smiles(row['smiles']),
    }))
    for _, row in df_D_X.iterrows()
}, namespace='D_X')
mapper.update({
    row['drug']: set(filter_none({
        row['drug'].strip().lower(),
        f"wikipedia:{row['drug'].strip().lower()}",
    }))
    for _, row in df_drugs_com.iterrows()
}, namespace='drugs.com')
mapper.update({
    'enalapril': {
        'enalapril',
        'enalaprilat',
    },
    'zofenopril-calcium': {
        'zofenopril-calcium',
        'zofenopril',
        'zofenopril calcium',
    },
}, namespace='manual')

mapper.update({
    f"CID:{int(row['CID'])}": set(filter_none({
        *extract_names(row['Drug (commercial)'].strip().lower() if type(row['Drug (commercial)']) == str else None),
        *extract_names(row['Preffered_name'].strip().lower() if not pd.isna(row['Preffered_name']) else None),
        f"drugcentral:{int(row['DrugCentral_ID'])}" if not pd.isna(row['DrugCentral_ID']) else None,
        rdkit_canonical_smiles(row['Smiles']) if row['Smiles'] else None,
    }))
    for _, row in df_placenta.iterrows()
    if not pd.isna(row['Cross']) and row['Cross'].lower().strip() == 'yes'
}, namespace='placenta')

mapper.summary()

#%%
df_placenta['name_mapped'] = [
    mapper.get(n, namespace='lincs_matrix')
    for i, n in enumerate(df_placenta['CID'].apply(lambda n: f"CID:{int(n)}"))
]
df_drugs_com['name_mapped'] =  [
    mapper.get(n, namespace='lincs_matrix')
    for i, n in enumerate(df_drugs_com['drug'].str.strip().str.lower())
]
df_D_X['name_mapped'] =  [
    mapper.get(n, namespace='lincs_matrix')
    for i, n in enumerate(df_D_X['name'].str.strip().str.lower())
]

# %%
s_P = set(df_placenta.loc[df_placenta['Cross'].str.lower()=='yes', 'name_mapped'].dropna())

s_X1 = set()
s_D1 = set()
for name, d in df_D_X.groupby('name_mapped'):
  if not name: continue
  if (d['category']=='X').any():
    s_X1.add(name)
  else:
    s_D1.add(name)

s_X2 = set()
s_D2 = set()
for name, d in df_drugs_com.groupby('name_mapped'):
  if not name: continue
  if (d['us-cleaned']=='X').any():
    s_X2.add(name)
  elif (d['us-cleaned']=='D').any():
    s_D2.add(name)

s_X = s_X1|s_X2
s_D = (s_D1|s_D2) - s_X

len(s_D), len(s_X), len(s_P)

# %%
supervenn(*zip(*[
  (s_X, 'Combined FDA Category X'), (s_X1, 'Drug Central FDA Category X'), (s_X2, 'Drugs.com FDA Category X'),
  (s_D, 'Combined FDA Category D'), (s_D1, 'Drug Central FDA Category D'), (s_D2, 'Drugs.com FDA Category D'),
]), widths_minmax_ratio=0.45,
    sets_ordering=None,
    chunks_ordering='minimize gaps',
)

# %%
supervenn(*zip(*[
  (set(coexpression_matrix_index), 'L1000 Perturbations'),
  (s_P, 'Placenta Crossing'),
  (s_D, 'FDA Category D'),
  (s_X, 'FDA Category X'),
]), widths_minmax_ratio=0.45,
    sets_ordering=None,
    chunks_ordering='minimize gaps',
)

plt.xlabel('DRUGS')
plt.savefig('data/2023-04-25-fig-2-supervenn.pdf')
plt.savefig('data/2023-04-25-fig-2-supervenn.png', dpi=300)

# %%
df_placenta.to_csv('data/2023-04-25-placenta-mapped.tsv', sep='\t')
df_D_X.to_csv('data/2023-04-25-D-X-mapped.tsv', sep='\t')
df_drugs_com.to_csv('data/2023-04-25-drugs-com.tsv', sep='\t')
