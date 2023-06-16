# %% [markdown]
# # ReproTox
# 
# - Data: FDA Pregnancy Categories D, X & Placental Crossing Set
# - Method: Visualize the drugs in the sigcom-lincs CD UMAP

# %%
import re
import pandas as pd
from supervenn import supervenn
from hdbscan import HDBSCAN
from matplotlib import pyplot as plt
from supervenn import supervenn
from maayanlab_bioinformatics.parse import gmt_read_dict

# %%
df_X_and_D = pd.read_excel('input/pregnancy_category_D_and_X_v2.xlsx')
df_X_and_D

# %%
df_placenta = pd.read_excel('input/Placenta_Barrier.xlsx')
df_placenta

# %%
df_X_and_D.groupby('name')['category'].nunique().value_counts()

# %%
s_X = set()
s_D = set()
for name, d in df_X_and_D.groupby('name'):
  if (d.category=='X').any():
    s_X.add(name)
  else:
    s_D.add(name)
len(s_D), len(s_X)

# %%
s_P = set(df_placenta[df_placenta['Cross'] == 'Yes']['Preffered_name'].str.strip().str.lower())

# %% [markdown]
# 

# %%
df_chem_meta = pd.read_csv('data/LINCS_small_molecules.tsv', sep='\t', index_col=0)
df_chem_meta

# %%
df_chem = pd.read_csv('data/2021-08-30-chem-umap.tsv.gz', sep='\t', index_col=0)
df_chem['cell'] = df_chem.index.map(lambda id: re.split(r'_', id)[1].split('.')[0])
df_chem['drug'] = df_chem.index.map(lambda id: re.split(r'_', id)[-2])
df_chem

# %%
s_l1000 = set(df_chem['drug'].str.lower()) & (s_D | s_X | s_P)
supervenn([s_D, s_X, s_P, s_l1000], ['Category D', 'Category X', 'Placenta', 'l1000'])

# %%
(s_D|s_X) - s_l1000

# %%
df_chem['FDA_category'] = df_chem['drug'].apply(lambda d: {
  (False, False): 'unknown',
  (False, True): 'Category D',
  (True, False): 'Category X',
}[(d in s_X, d in s_D)])

# %%
df_chem['Placenta'] = df_chem['drug'].apply(lambda d: {False: 'unknown', True: 'Placenta Crossing'}[d in s_P])

# %%
df_chem['drug'].isin(list(s_l1000)).value_counts()

# %%
moa_to_drug = gmt_read_dict(
    'data/DrugRepurposingHub_moa_drugsetlibrary_name.dmt',
    parse_gene=lambda r: (r, 1)
)
drug_to_moa = pd.Series({
    drug: moa
    for moa, drugs in moa_to_drug.items()
    for drug in drugs
}).to_frame('MOA')
drug_to_moa

# %%
df_chem_w_moa = pd.merge(left=df_chem, left_on='drug', right=drug_to_moa, right_index=True).fillna('unknown')
df_chem_w_moa

# %%
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors

# %%
fig,ax=plt.subplots(figsize=(14,12))
D = {k: d for k, d in df_chem.groupby('FDA_category')}
cat = 'unknown'
d = D.pop(cat)
ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey', rasterized=True)
for cat, d in D.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=1, alpha=0.5, label=cat, rasterized=True)
plt.xlim((3, 16))
plt.ylim((3, 16))
# l = plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
l = plt.legend(loc='lower right', ncol=2, fontsize=14)
# l = plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize=14)
for h in l.legendHandles: h._sizes = [50]
plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plt.tight_layout()
plt.savefig('figures/l1000-space-fda.svg', bbox_inches='tight', dpi=300)
plt.savefig('figures/l1000-space-fda.png', bbox_inches='tight', dpi=300)
plt.show()

# %%
fig,ax=plt.subplots(figsize=(14,12))
D = {k: d for k, d in df_chem.groupby('Placenta')}
cat = 'unknown'
d = D.pop(cat)
ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey', rasterized=True)
for cat, d in D.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=1, alpha=0.5, label=cat, c='#2ca02c', rasterized=True)
plt.xlim((3, 16))
plt.ylim((3, 16))
# l = plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
l = plt.legend(loc='lower right', ncol=2, fontsize=14)
# l = plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize=14)
for h in l.legendHandles: h._sizes = [50]
plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plt.tight_layout()
plt.savefig('data/fig-3-l1000-space-placenta.pdf', bbox_inches='tight', dpi=300)
plt.savefig('data/fig-3-l1000-space-placenta.png', bbox_inches='tight', dpi=300)
plt.show()

# %%
fig,ax=plt.subplots(figsize=(12,6))
d = df_chem_w_moa[((df_chem_w_moa['UMAP-1'] > 9) & (df_chem_w_moa['UMAP-1'] < 12) & (df_chem_w_moa['UMAP-2'] > 5) & (df_chem_w_moa['UMAP-2'] < 7))]
D = {k: d for k, d in d.groupby('MOA') if d.shape[0]>800}
D_other = {k: d for k, d in d.groupby('MOA') if d.shape[0]<=800}
# cat = 'unknown'
for cat, d in D_other.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=1, c='lightgrey')
for cat, d in D.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=10, alpha=0.5, label=cat)
# plt.xlim((9, 12))
# plt.ylim((5, 7))
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.tight_layout()
plt.show()

# %%
fda_drug_moas = df_chem_w_moa.groupby('FDA_category')['MOA'].value_counts().unstack().drop('unknown').dropna(how='all', axis=1).T.sort_values('Category D')
fda_drug_moas

# %%
fda_top_drug_moas = fda_drug_moas[fda_drug_moas.rank(ascending=False) < 5].dropna(how='all')
fda_top_drug_moas

# %%
fig,ax=plt.subplots(figsize=(14,12))
D = {k: d for k, d in df_chem_w_moa.groupby('MOA') if k in fda_top_drug_moas.index}
D_other = {k: d for k, d in df_chem_w_moa.groupby('MOA') if k not in fda_top_drug_moas.index}
for cat, d in D_other.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey')
# cat = 'unknown'
# d = D.pop(cat)
# ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey')
for cat, d in D.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=2, alpha=0.5, label=cat)
plt.xlim((3, 16))
plt.ylim((3, 16))
l = plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
for h in l.legendHandles: h._sizes = [30]
plt.tight_layout()
# plt.savefig('figures/pregnancy-map-moa.png', dpi=300)
plt.show()

# %%
clusterer = HDBSCAN(min_cluster_size=40)
df_chem_w_moa['clust'] = clusterer.fit_predict(df_chem_w_moa.loc[:, ['UMAP-1', 'UMAP-2']])

# %%
df_chem_w_moa['clust'].value_counts()

# %%
fig,ax=plt.subplots(figsize=(14,12))
D = {k: d for k, d in df_chem_w_moa.groupby('clust')}
# cat = 'unknown'
# d = D.pop(cat)
# ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey')
for cat, d in D.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=2, alpha=0.5, label=cat)
plt.xlim((3, 16))
plt.ylim((3, 16))
# l = plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
# for h in l.legendHandles: h._sizes = [30]
plt.tight_layout()
plt.show()

# %%
df_chem_w_moa

# %%
d = pd.get_dummies(df_chem_w_moa['clust']).corrwith(pd.get_dummies(df_chem_w_moa['FDA_category']).loc[:, 'Category X']).rank(ascending=False)
df_chem_w_moa['clust_X'] = df_chem_w_moa['clust'].isin(d[d < 25].index)
d = pd.get_dummies(df_chem_w_moa['clust']).corrwith(pd.get_dummies(df_chem_w_moa['FDA_category']).loc[:, 'Category D']).rank(ascending=False)
df_chem_w_moa['clust_D'] = df_chem_w_moa['clust'].isin(d[d < 25].index)
d = pd.get_dummies(df_chem_w_moa['clust']).corrwith(pd.get_dummies(df_chem_w_moa['Placenta']).loc[:, 'Placenta Crossing']).rank(ascending=False)
df_chem_w_moa['clust_P'] = df_chem_w_moa['clust'].isin(d[d < 25].index)

# %%
c = df_chem_w_moa.loc[df_chem_w_moa['clust_X'], 'MOA'].value_counts()
clust_X_top_moas = c.iloc[:5]
c = df_chem_w_moa.loc[df_chem_w_moa['clust_D'], 'MOA'].value_counts()
clust_D_top_moas = c.iloc[:5]
c = df_chem_w_moa.loc[df_chem_w_moa['clust_P'], 'MOA'].value_counts()
clust_P_top_moas = c.iloc[:5]

# %%
supervenn([set(fda_top_drug_moas.index), set(clust_X_top_moas.index), set(clust_D_top_moas.index), set(clust_P_top_moas.index)], ['fda', 'X', 'D', 'P'])

# %%
set(fda_top_drug_moas.index) - set(clust_D_top_moas.index) - set(clust_X_top_moas.index)

# %%
(set(clust_D_top_moas.index) | set(clust_X_top_moas.index)) - set(fda_top_drug_moas.index)

# %%
fig,ax=plt.subplots(figsize=(14,12))
D = {k: d for k, d in df_chem_w_moa.groupby('MOA') if k in set(clust_D_top_moas.index) | set(clust_X_top_moas.index)}# | set(clust_P_top_moas.index)}
D_other = {k: d for k, d in df_chem_w_moa.groupby('MOA') if k not in set(clust_D_top_moas.index) | set(clust_X_top_moas.index)}# | set(clust_P_top_moas.index)}
d = df_chem[~df_chem.index.isin(df_chem_w_moa)]
ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey', rasterized=True)
for cat, d in D_other.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey', rasterized=True)
# cat = 'unknown'
# d = D.pop(cat)
# ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=0.1, c='lightgrey')
for cat, d in D.items():
  ax.scatter(x=d['UMAP-1'], y=d['UMAP-2'], s=2, alpha=0.5, label=cat.replace(' (SERM)', ''), rasterized=True)
plt.xlim((3, 16))
plt.ylim((3, 16))
# l = plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
l = plt.legend(loc='lower right', ncol=2, fontsize=14)
for h in l.legendHandles: h._sizes = [50]
plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
plt.tight_layout()
plt.savefig('figures/fig-3-l1000-space-moa.svg', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig-3-l1000-space-moa.png', dpi=300, bbox_inches='tight')
plt.show()


