# %%
import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.metrics import pairwise_distances

# %%
rdkit_features = pd.HDFStore('data/drug_attribute_mats.h5')

# %%
all_features = rdkit_features['all_features']
all_features

# %%
groups = {}
for index, row in all_features.iterrows():
  features = row.to_numpy()
  b = features.data.tobytes()
  if b not in groups:
    groups[b] = set()
  groups[b].add(index)

# %%
groups_sorted = sorted([(len(v), k) for k, v in groups.items()], reverse=True)

# %%
all_features_reduced = all_features.drop_duplicates().droplevel(0)
all_features_reduced

# %%
rdkit_features['all_features_reduced'] = all_features_reduced

# %%
drug_to_smile = {drug: smile for V in groups.values() for smile, drug in V}
in_reduced = {
  (drug_to_smile[drug], drug)
  for drug in all_features_reduced.index
}
drug_syn = {
  next(iter(V & in_reduced))[1]: {drug for smile, drug in V}
  for V in groups.values()
}

#%%
all_features_reduced_idf = all_features_reduced/ all_features_reduced.sum()
all_features_reduced_idf

# %%
all_features_reduced_idf_sim = pd.DataFrame(
  1. - pairwise_distances(all_features_reduced_idf.values, metric='cosine', n_jobs=4),
  index=all_features_reduced_idf.index, columns=all_features_reduced_idf.index)
np.fill_diagonal(all_features_reduced_idf_sim.values, 0)

# %%
rdkit_features['all_features_reduced_idf_sim'] = all_features_reduced_idf_sim

# %%
rdkit_feature_types = all_features_reduced_idf.columns.map(lambda s: s.split('_')[0]).value_counts()
cols = all_features_reduced_idf.columns
for t in rdkit_feature_types.index:
  if f'{t}_reduced_idf_sim' in rdkit_features: continue
  print(cols[cols.str.startswith(t + '_')])
  sim = pd.DataFrame(
    1. - pairwise_distances(
      all_features_reduced_idf.loc[:, cols[cols.str.startswith(t + '_')]].values,
      metric='cosine',
      n_jobs=4,
    ),
    index=all_features_reduced_idf.index,
    columns=all_features_reduced_idf.index,
  )
  np.fill_diagonal(sim.values, 0)
  rdkit_features[f'{t}_reduced_idf_sim'] = sim

#%%
physio_props = rdkit_features['physio_props']
physio_props_reduced = physio_props.drop_duplicates().droplevel(0)
physio_props_reduced_zscore = zscore(physio_props_reduced)
physio_props_reduced_zscore_sim = pd.DataFrame(
  1. - pairwise_distances(
    physio_props_reduced_zscore,
    metric='correlation',
    n_jobs=4,
  ),
  index=physio_props_reduced_zscore.index,
  columns=physio_props_reduced_zscore.index,
)
np.fill_diagonal(physio_props_reduced_zscore_sim.values, 0)
rdkit_features['physio_props_reduced_zscore_sim'] = physio_props_reduced_zscore_sim

#%%
for t in physio_props_reduced_zscore.columns:
  k = f'{t}_physio_props_reduced_zscore_sim'
  # if k in rdkit_features: continue
  sim = pd.DataFrame(
    1. - pairwise_distances(
      physio_props_reduced_zscore.loc[:, [t]].values,
      metric='correlation',
      n_jobs=4,
    ),
    index=physio_props_reduced_zscore.index,
    columns=physio_props_reduced_zscore.index,
  )
  np.fill_diagonal(sim.values, 0)
  rdkit_features[k] = sim

# %%
rdkit_features.close()
