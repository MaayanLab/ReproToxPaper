# %%
import json
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import MACCSkeys
from pathlib import Path

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

def rdkit_canonical_smiles(smiles):
    if not smiles or smiles == '-': return None
    try:
        RDLogger.DisableLog('rdApp.*')
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    except KeyboardInterrupt:
        raise
    except:
        return smiles


# %%
import uuid
from collections import Counter

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
            self._forward[id] = Counter()
            for synonym in (key, *synonyms):
                if synonym not in self._reverse:
                    self._forward[id].update([synonym])
                    self._reverse[synonym] = id
                    if synonym not in self._namespaces:
                        self._namespaces[synonym] = {}
                    if namespace not in self._namespaces[synonym]:
                        self._namespaces[synonym][namespace] = Counter()
                    self._namespaces[synonym][namespace].update([key])
                else:
                    orig_id = self._reverse[synonym]
                    if orig_id != id:
                        for s in self._forward[id]:
                            self._forward[orig_id].update([s])
                            self._reverse[s] = orig_id
                        if synonym not in self._namespaces:
                            self._namespaces[synonym] = {}
                        if namespace not in self._namespaces[synonym]:
                            self._namespaces[synonym][namespace] = Counter()
                        self._namespaces[synonym][namespace].update([key])
                        del self._forward[id]
                        id = orig_id


# %%
coexpression_matrix = np.load('data/L1000_2021_drug_similarity.npz', allow_pickle=True)
coexpression_matrix_index = coexpression_matrix['index'][:]

# %%
lincs_sm_meta = pd.read_csv('data/LINCS_small_molecules.tsv', sep='\t', index_col=0)
lincs_sm_meta

# %%
drugbank_synonyms = read_gmt_dict('data/drugbank_drug_synonyms.dmt')
drugbank_enzyme = pd.DataFrame({k: {v: 1 for v in V} for k, V in read_gmt_iter('data/drugbank_enzyme.dmt')})
drugbank_transporter = pd.DataFrame({k: {v: 1 for v in V} for k, V in read_gmt_iter('data/drugbank_transporter.dmt')})
drugshot_autorif_enzyme_drugbank = pd.read_csv('data/drugshot_autorif_enzyme_drugbank.tsv', sep='\t', index_col=0)
drugshot_autorif_transporter_drugbank = pd.read_csv('data/drugshot_autorif_transporter_drugbank.tsv', sep='\t', index_col=0)
drugshot_drugrif_enzyme_drugbank = pd.read_csv('data/drugshot_drugrif_enzyme_drugbank.tsv', sep='\t', index_col=0)
drugshot_drugrif_transporter_drugbank = pd.read_csv('data/drugshot_drugrif_transporter_drugbank.tsv', sep='\t', index_col=0)
chembl_drugs = {
    drug['chemblId']: drug
    for drug in json.load(Path('data/chembl_drugs.json').open('r'))
}

# %%
mapper = IDMapper()
mapper.update({
    i: {
        s.strip().lower()
        for s in [
            i,
            row['pert_name'],
            row['canonical_smiles'],
            row['inchi_key'],
        ]
        if s
    }
    for i, row in lincs_sm_meta.iterrows()
}, namespace='lincs_sm_meta')
mapper.update({
  drug['chemblId']: {
      synonym.strip().lower()
      for synonym in (drug['synonyms'] + [
          drug['prefName'],
          f"ChEMBL:{drug['chemblId']}",
      ])
      if synonym
  }
  for drug in chembl_drugs.values()
}, namespace='chembl')
mapper.update({
    d: {
        s.strip().lower()
        for s in S
        if s
    }
    for d, S in drugbank_synonyms.items()
}, namespace='drugbank')
mapper.update({
    rdkit_canonical_smiles(row['canonical_smiles']): {i}
    for i, row in lincs_sm_meta.iterrows()
}, namespace='smiles')
mapper.update({
    rdkit_canonical_smiles(drug['smiles']): {i}
    for i, drug in chembl_drugs.items()
}, namespace='smiles')
mapper.update({
    idx: {idx.strip().lower()}
    for idx in coexpression_matrix_index
}, namespace='lincs_matrix')
mapper.summary()

# %%
all_drugs = coexpression_matrix_index
all_drug_smiles = {drug: smiles for drug in all_drugs for smiles_ in (mapper.get(drug, namespace='smiles'),) if smiles_ for smiles in smiles_.keys()}
all_mols = {smiles: Chem.MolFromSmiles(smiles) for smiles in set(all_drug_smiles.values())}

# %%
all_smiles_drugs = {}
for drug, smiles in all_drug_smiles.items():
    if smiles not in all_smiles_drugs:
        all_smiles_drugs[smiles] = set()
    all_smiles_drugs[smiles].add(drug)

# %%
from tqdm import tqdm
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

# %%
def extract(prefix, v):
    return {f"{prefix}_{i}": np.int8(1) for i in v.GetOnBits()}

all_features = []
for smiles, mol in tqdm(all_mols.items(), total=len(all_mols)):
    info = dict(
        smiles=smiles,
        **extract('maccs', MACCSkeys.GenMACCSKeys(mol)),
        **extract('morgan2f', AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048, useFeatures=True)),
        **extract('morgan3f', AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048, useFeatures=True)),
        **extract('morgan4f', AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=2048, useFeatures=True)),
        **extract('avalon_', pyAvalonTools.GetAvalonFP(mol)),
        **extract('atomPair', Chem.rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)),
        **extract('rdkFingerprint2', Chem.RDKFingerprint(mol, maxPath=2)),
        **extract('rdkFingerprint5', Chem.RDKFingerprint(mol, maxPath=5)),
        **extract('topological', FingerprintMols.FingerprintMol(mol, minPath=1, maxPath=7, fpSize=2048, bitsPerHash=2, useHs=True, tgtDensity=0, minSize=128)),
    )
    for drug in all_smiles_drugs[smiles]:
        all_features.append(dict(drug=drug, **info))

# %%
df_all_features = pd.DataFrame(all_features).set_index(['smiles', 'drug']).fillna(0).astype(np.int8)
df_all_features

# %%
store = pd.HDFStore('data/drug_attribute_mats.h5')
store['all_features'] = df_all_features
store.close()

# %% [markdown]
# # Physio-Chemical Properties

# %%
from rdkit.Chem import QED
from rdkit.Chem import Crippen

# %%
physio_props = []
for smiles, mol in tqdm(all_mols.items(), total=len(all_mols)):
    info = dict(
        smiles=smiles,
        **{f"qed_{prop}": value for prop, value in QED.properties(mol)._asdict().items()},
        crippen_logp=Crippen.MolLogP(mol),
        crippen_mr=Crippen.MolMR(mol),
    )
    for drug in all_smiles_drugs[smiles]:
        physio_props.append(dict(drug=drug, **info))

# %%
df_physio_props = pd.DataFrame(physio_props).set_index(['smiles', 'drug'])
df_physio_props

# %%
store = pd.HDFStore('data/drug_attribute_mats.h5')
store['physio_props'] = df_physio_props
store.close()
