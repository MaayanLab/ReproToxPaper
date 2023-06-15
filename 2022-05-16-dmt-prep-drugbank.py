# %% [markdown]
# # CFDE Toxicology Partnership: Drug Matrix Transpose (DMT) Preparation
# 
# Prepare enzyme/transporter => drug associations as DMTs from various sources.
# 
# - [x] Metrabase
# - [x] ChEMBL<sup>1</sup>
# - [x] DrugBank
# - [x] DrugShot PubMed Queries (AutoRif & DrugRif)
# 
# 
# <sup>1</sup> only metabolism currently

# %%
import time
import pandas as pd
import requests
import traceback
import json
from pathlib import Path

# %%
def call_with_exp_backoff(func, *args, **kwargs):
    backoff = 1
    while True:
        try:
            return func(*args, **kwargs)
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()
        time.sleep(backoff)
        backoff *= 2

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

def ensure_list(L):
    return L if type(L) == list else [L]

def one(it):
    try:
        return next(iter(it))
    except TypeError:
        return None
    except StopIteration:
        return None

def read_gmt_iter(p: Path):
    with Path(p).open('r') as fr:
        for line in fr:
            term, *geneset = line.strip('\n').split('\t')
            yield term, set(filter(None, geneset))

def read_gmt_dict(p: Path):
    return dict(read_gmt_iter(p))


def one(it):
    try:
        return next(iter(it))
    except TypeError:
        return None
    except StopIteration:
        return None

import uuid

class IDMapper:
    ''' Stores id mappings and makes it easy to use many of them in tandem
    '''
    def __init__(self):
        self._forward = {}
        self._reverse = {}
        self._namespaces = {}
        self._counts = {}

    def summary(self):
        summary = {}
        for synonyms in self._forward.values():
            ns = frozenset(
                ns
                for synonym in synonyms
                for ns in self._namespaces[synonym]
            )
            if ns not in summary:
                summary[ns] = 0
            summary[ns] += 1
        return summary
        
    def get(self, term, namespace=None):
        id = self._reverse.get(term)
        if id is None: return None
        refs = {
            namespace: keys
            for synonym in self._forward[id]
            for namespace, keys in self._namespaces[synonym].items()
        }
        return dict(
            id=id,
            refs=refs,
            synonyms=self._forward[id],
        ) if namespace is None else refs.get(namespace)

    def update(self, mappings, namespace=None):
        ''' Add mappings of the form:
        { identifier: { synonyms } }
        '''
        for key, synonyms in mappings.items():
            id = uuid.uuid4()
            self._forward[id] = set()
            for synonym in (key, *synonyms):
                if synonym not in self._reverse:
                    self._forward[id].add(synonym)
                    self._reverse[synonym] = id
                    if synonym not in self._namespaces:
                        self._namespaces[synonym] = {}
                    if namespace not in self._namespaces[synonym]:
                        self._namespaces[synonym][namespace] = set()
                    self._namespaces[synonym][namespace].add(key)
                else:
                    orig_id = self._reverse[synonym]
                    if orig_id != id:
                        for s in self._forward[id]:
                            self._forward[orig_id].add(s)
                            self._reverse[s] = orig_id
                        if synonym not in self._namespaces:
                            self._namespaces[synonym] = {}
                        if namespace not in self._namespaces[synonym]:
                            self._namespaces[synonym][namespace] = set()
                        self._namespaces[synonym][namespace].add(key)
                        del self._forward[id]
                        id = orig_id


# %% [markdown]
# ## DrugBank

# %%
d = Path('data')
drugbank_synonyms = read_gmt_dict(d/'drugbank_drug_synonyms.dmt')
drugbank_enzyme = read_gmt_dict(d/'drugbank_enzyme.dmt')
drugbank_transporter = read_gmt_dict(d/'drugbank_transporter.dmt')

# %% [markdown]
# ## DrugShot Pubmed Search
# 
# ![image.png](attachment:51a73295-c4dd-4c64-a703-b52f07ccdaae.png)

# %%
class NoResults(Exception): pass
def drugshot_associate(term, rif='autorif'):
  ret = call_with_exp_backoff(requests.post, 'https://maayanlab.cloud/drugshot/api/search', json=dict(rif=rif, term=term)).json()
  if not ret['drug_count']: raise NoResults
  return pd.DataFrame({ drug: dict(count=count, frac=frac, score=count * frac) for drug, [count,frac] in ret['drug_count'].items() }).T.sort_values('score', ascending=False)

# %%
import re
expr = re.compile(r'^(.+) (.+?)$')

# %%
enzyme_terms = {expr.match(e).group(1) for e, ds in drugbank_enzyme.items() if len(ds) >= 5}
transporter_terms = {expr.match(e).group(1) for e, ds in drugbank_transporter.items() if len(ds) >= 5}

# %%
drugshot_autorif_enzyme = {}
for symbol in enzyme_terms:
  if symbol in drugshot_autorif_enzyme: continue
  try:
    drugshot_autorif_enzyme[symbol] = drugshot_associate(symbol, rif='autorif').to_dict()
  except NoResults:
    pass

# %%
with open('data/drugshot_autorif_enzyme.json', 'w') as fw:
    json.dump(drugshot_autorif_enzyme, fw)

# %%
drugshot_autorif_transporter = {}
for symbol in transporter_terms:
  if symbol in drugshot_autorif_transporter: continue
  try:
    drugshot_autorif_transporter[symbol] = drugshot_associate(symbol, rif='autorif').to_dict()
  except NoResults:
    pass

# %%
with open('data/drugshot_autorif_transporter.json', 'w') as fw:
    json.dump(drugshot_autorif_transporter, fw)

# %%
drugshot_drugrif_enzyme = {}
for symbol in enzyme_terms:
  if symbol in drugshot_drugrif_enzyme: continue
  try:
    drugshot_drugrif_enzyme[symbol] = drugshot_associate(symbol, rif='drugrif').to_dict()
  except NoResults:
    pass

# %%
with open('data/drugshot_drugrif_enzyme.json', 'w') as fw:
    json.dump(drugshot_drugrif_enzyme, fw)

# %%
drugshot_drugrif_transporter = {}
for symbol in transporter_terms:
  if symbol in drugshot_drugrif_transporter: continue
  try:
    drugshot_drugrif_transporter[symbol] = drugshot_associate(symbol, rif='drugrif').to_dict()
  except NoResults:
    pass

# %%
with open('data/drugshot_drugrif_transporter.json', 'w') as fw:
    json.dump(drugshot_drugrif_transporter, fw)

# %%


# %%
pd.DataFrame({k: v['score'] for k,v in drugshot_drugrif_transporter.items()}).to_csv('data/drugshot_drugrif_transporter_drugbank.tsv', sep='\t')

# %%
pd.DataFrame({k: v['score'] for k,v in drugshot_autorif_transporter.items()}).to_csv('data/drugshot_autorif_transporter_drugbank.tsv', sep='\t')

# %%
pd.DataFrame({k: v['score'] for k,v in drugshot_drugrif_enzyme.items()}).to_csv('data/drugshot_drugrif_enzyme_drugbank.tsv', sep='\t')

# %%
pd.DataFrame({k: v['score'] for k,v in drugshot_autorif_enzyme.items()}).to_csv('data/drugshot_autorif_enzyme_drugbank.tsv', sep='\t')


