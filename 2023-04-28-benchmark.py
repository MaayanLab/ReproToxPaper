# %%
import glasbey
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow_addons as tfa
import multiprocessing as mp
import scipy.stats
from tqdm.auto import tqdm

def nan_diag(m, copy=False):
    if copy: m = m.copy()
    np.fill_diagonal(m.values, np.nan)
    return m

# %%
df_placenta = pd.read_csv('data/2023-04-25-placenta-mapped.tsv', sep='\t')
df_D_X = pd.read_csv('data/2023-04-25-D-X-mapped.tsv', sep='\t')
df_drugs_com = pd.read_csv('data/2023-04-25-drugs-com.tsv', sep='\t')

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
rdkit_features = pd.HDFStore('data/drug_attribute_mats.h5', 'r')
all_features_reduced = rdkit_features['all_features_reduced']
rdkit_feature_types = all_features_reduced.columns.map(lambda s: s.split('_')[0]).value_counts()
del all_features_reduced

# %%
coexpression_matrix = np.load('data/L1000_2021_drug_similarity.npz', allow_pickle=True)
coexpression_matrix_index = coexpression_matrix['index'][:]
coexpression_matrix = pd.DataFrame(
    data=coexpression_matrix['correlations'],
    columns=coexpression_matrix_index,
    index=coexpression_matrix_index,
)
nan_diag(coexpression_matrix)

# %%
backgrounds = {
  'All Structure Similarity': lambda: rdkit_features[f"all_features_reduced_idf_sim"],
  'L1000 Similarity': lambda: coexpression_matrix,
  'Rescale[L1000 Similarity]': lambda: (coexpression_matrix + 1) / 2,
  'physio_props': lambda: rdkit_features['physio_props_reduced_zscore_sim'],
  'Rescale[physio_props]': lambda: (rdkit_features['physio_props_reduced_zscore_sim'] + 1)/2,
}
backgrounds.update({
  t: lambda t=t: rdkit_features[f"{t}_reduced_idf_sim"]
  for t in rdkit_feature_types.index
})

# %%
def train_test_split(L, random_seed=42):
    L = list(L)
    L.sort()
    np.random.seed(random_seed)
    np.random.shuffle(L)
    return L[len(L)//3:], L[:len(L)//3]

labels = {
  'FDA Category D': train_test_split(s_D),
  'FDA Category X': train_test_split(s_X),
  'FDA Category D & X': train_test_split(s_D|s_X),
  'Placenta Crossing': train_test_split(s_P),
}

# %%
def shuffle(x):
    x = x.copy()
    np.random.shuffle(x)
    return x

def enrichment_score(y_score, y_true):
    hits, = np.where(np.in1d(y_score.index, y_true))
    hit_indicator = np.zeros(y_score.shape[0])
    hit_indicator[hits] = 1
    no_hit_indicator = 1 - hit_indicator
    number_hits = len(hits)
    number_miss = y_score.shape[0] - number_hits
    sum_hit_scores = np.sum(y_score.iloc[hits].abs())
    if sum_hit_scores == 0:
        return 0, hits, 0
    norm_hit =  float(1.0/sum_hit_scores)
    norm_no_hit = float(1.0/number_miss)
    running_sum = np.cumsum(hit_indicator * y_score.abs().values * norm_hit - no_hit_indicator * norm_no_hit)
    nn = np.argmax(np.abs(running_sum))
    es = running_sum[nn]
    return running_sum, hits, es

def enrichment_score_star(args):
    return enrichment_score(*args)

def normalize_enrichment_score(y_score, y_true, permutations=10000):
    running_sum, hits, es = enrichment_score(y_score, y_true)
    running_sum_nulls = []
    esnulls = []
    for running_sum_null, _, esnull in tqdm(map(
        enrichment_score_star,
        (
            (pd.Series(y_score.values, index=shuffle(y_score.index.values)), y_true)
            for _ in range(permutations)
        )
    ), total=permutations):
        running_sum_nulls.append(running_sum_null)
        esnulls.append(esnull)
    running_sum_nulls = np.array(running_sum_nulls)
    running_sum_null_mu = running_sum_nulls.mean(axis=0)
    running_sum_null_std = running_sum_nulls.std(axis=0)
    esnulls = np.array(esnulls)
    if es < 0:
        pval = (esnulls < es).sum() / (esnulls < 0).sum()
        esnull_mu = esnulls[esnulls<0].mean()
        esnull_std = esnulls[esnulls<0].std()
        nes = -es / esnull_mu
    elif es >= 0:
        pval = (esnulls >= es).sum() / (esnulls >= 0).sum()
        esnull_mu = esnulls[esnulls>=0].mean()
        esnull_std = esnulls[esnulls>=0].std()
        nes = es / esnull_mu
    z_es = (es-esnull_mu)/esnull_std
    z_es_pval = scipy.stats.norm.sf(abs(z_es))*2
    return dict(
        running_sum=running_sum,
        hits=hits,
        es=es,
        esnull_mu=esnull_mu,
        esnull_std=esnull_std,
        nes=nes,
        pval=pval,
        z_es=z_es,
        z_es_pval=z_es_pval,
        running_sum_null_mu=running_sum_null_mu,
        running_sum_null_std=running_sum_null_std,
    )

# %%
def flatten(y_score):
    ''' Uniform weights, ordering preserved
    '''
    return pd.Series(np.ones(y_score.shape[0]), index=y_score.index)

def get_scores_for(background_matrix, drugs):
    background_matrix = background_matrix()
    mat = background_matrix.loc[np.in1d(background_matrix.index, list(drugs))].copy()
    del background_matrix
    for drug in mat.index: mat.loc[drug, drug] = np.nan
    mat.loc['Similarity_Score'] = mat[mat.columns].mean()
    return mat.loc['Similarity_Score'].sort_values(ascending = False)

def benchmark(args):
    label, (train_drugs, test_drugs), background = args
    y_score = get_scores_for(backgrounds[background], train_drugs)
    y_score_test = y_score[~y_score.index.isin(train_drugs)]
    return dict(
        normalize_enrichment_score(flatten(y_score_test), test_drugs),
        y_score=y_score,
        label=label,
        background=background,
    )

# %%
with mp.Pool() as p:
    results = pd.DataFrame(list(tqdm(
        p.imap_unordered(benchmark, (
            (label, train_test_drugs, background)
            for label, train_test_drugs in labels.items()
            for background in backgrounds
        )),
        total=len(backgrounds)*len(labels),
    )))

# %%
def linear_combine_score_model(X, train_drugs):
    y_train = pd.Series(X.index.isin(train_drugs).astype(float), index=X.index)
    r = y_train.sum() / y_train.shape[0]
    m = tf.keras.models.Sequential([
        tf.keras.layers.Input((X.shape[1],)),
        tf.keras.layers.Dense(1, activation='sigmoid'),
    ])
    m.compile(optimizer='adam', loss='binary_crossentropy')
    m.fit(X.values, y_train.values,
          epochs=50, class_weight={ 0: 1-r, 1: r },
          verbose=0,
          callbacks=[
            tfa.callbacks.TQDMProgressBar(show_epoch_progress=False),
          ])
    y_score = pd.Series(np.squeeze(m(X.values)), index=X.index)
    # lr = LogisticRegressionCV(class_weight='balanced', random_state=42, max_iter=1000)
    # lr.fit(X, y_train)
    # y_score = pd.Series(lr.predict_proba(X)[:, 1], index=X.index)
    return y_score

def benchmark_agg(args):
    label, (train_drugs, test_drugs), background, y_score = args
    y_score_test = y_score[~y_score.index.isin(train_drugs)]
    return dict(
        normalize_enrichment_score(flatten(y_score_test), test_drugs),
        y_score=y_score,
        label=label,
        background=background,
    )

def reweight(y_score):
    r = pd.Series(y_score).rank(ascending=True, method='dense') - 1
    return (r - r.min()) / (r.max() - r.min())

tasks = []
for bgs in tqdm([
    ['All Structure Similarity', 'physio_props', 'L1000 Similarity'],
    ['physio_props', 'L1000 Similarity'],
    ['morgan4f', 'L1000 Similarity'],
    ['morgan4f', 'physio_props'],
    ['morgan4f', 'physio_props', 'L1000 Similarity'],
]):
    # pre_combined_model = precompute_model(bgs)
    for label, (train_drugs, test_drugs) in tqdm(labels.items()):
        records = results.loc[((results.label==label))]
        x = pd.concat(records.loc[records.background.isin(bgs), 'y_score'].tolist(), axis=1)
        top_rank = 1-reweight(x.rank(ascending=False).min(axis=1))
        mean_rank = 1-reweight(x.rank(ascending=False).mean(axis=1))
        x = x.dropna()
        x = (x - x.mean()) / np.maximum(x.std(), 1e-9)
        linear_combined_score = linear_combine_score_model(x, train_drugs)
        tasks.append((label, (train_drugs, test_drugs), f"top_rank[{','.join(set(bgs))}]", top_rank))
        tasks.append((label, (train_drugs, test_drugs), f"mean_rank[{','.join(set(bgs))}]", mean_rank))
        tasks.append((label, (train_drugs, test_drugs), f"linear_combined[{','.join(set(bgs))}]", linear_combined_score))

with mp.Pool() as p:
    agg_results = pd.DataFrame(list(tqdm(
        p.imap_unordered(benchmark_agg, tasks),
        total=len(tasks),
    )))

# %%
C = ['label', 'background', 'es', 'nes', 'pval', 'z_es', 'z_es_pval', 'esnull_mu', 'esnull_std', 'running_sum', 'running_sum_null_mu', 'running_sum_null_std']
all_results = pd.concat([results[C], agg_results[C]])

# %%
all_results['ks'] = all_results.apply(lambda row: scipy.stats.ks_2samp(row['running_sum'], row['running_sum_null_mu']), axis=1)
all_results['ks_stat'] = all_results['ks'].apply(lambda ks: ks.statistic)
all_results['ks_p'] = all_results['ks'].apply(lambda ks: ks.pvalue)

# %%
all_results.sort_values('ks_stat', ascending=False)[['label', 'background', 'nes', 'z_es', 'ks_stat', 'ks_p']].iloc[:40]

# %%
# row = all_results.loc[((all_results.label=='FDA Category D & X')&(all_results.background=='rdkFingerprint2')&(~all_results.reweighted))].iloc[0]
row = all_results.loc[((all_results.label=='Placenta Crossing')&(all_results.background=='Rescale[L1000 Similarity]'))].iloc[0]
plt.plot(
    np.arange(row['running_sum'].shape[0]) / row['running_sum'].shape[0],
    row['running_sum'],
)
plt.plot(
    np.arange(row['running_sum_null_mu'].shape[0]) / row['running_sum_null_mu'].shape[0],
    row['running_sum_null_mu'] + 2*row['running_sum_null_std'],
    linestyle='dashdot',
)
plt.plot(
    np.arange(row['running_sum_null_mu'].shape[0]) / row['running_sum_null_mu'].shape[0],
    row['running_sum_null_mu'],
    linestyle='dashed',
)
plt.plot(
    np.arange(row['running_sum_null_mu'].shape[0]) / row['running_sum_null_mu'].shape[0],
    row['running_sum_null_mu'] - 2*row['running_sum_null_std'],
    linestyle='dashdot',
)

# %%
#%%
import re
# bgs = all_results.sort_values('z_es')['background'].iloc[-20:].unique().tolist() + ['Rescale[L1000 Similarity]', 'Rescale[physio_props]', 'morgan4f']
bgs = [
    'morgan4f',
    'physio_props',
    'L1000 Similarity',
    'mean_rank[morgan4f,physio_props,L1000 Similarity]',
    'top_rank[morgan4f,physio_props,L1000 Similarity]',
    'linear_combined[morgan4f,physio_props,L1000 Similarity]',
    # 'deep_combined[physio_props,morgan4f,L1000 Similarity]',
]
colors = dict(zip(bgs, glasbey.create_palette(len(bgs), lightness_bounds=(0, 64))))

#%%
#%%
for label, records in all_results.groupby('label'):
    if label not in {'FDA Category D & X', 'Placenta Crossing'}: continue
    fig, ((ax11, ax12, ax13, ax14), (ax21, ax22, ax23, ax24), (ax31, ax32, ax33, ax34), (ax41, ax42, ax43, ax44)) = plt.subplots(4, 4, figsize=(25, 25))
    for _, row in records.iterrows():
        if row['background'] not in colors: continue
        running_sum = row['running_sum'].copy()
        running_sum_null_mu = row['running_sum_null_mu'].copy()
        running_sum_null_mu[-1] = float('nan')
        running_sum_null_std = row['running_sum_null_std'].copy()
        running_sum_null_std[-1] = float('nan')
        running_sum_z = (running_sum - running_sum_null_mu) / np.maximum(running_sum_null_std, 1e-9)
        running_sum_z[-1] = float('nan')
        running_sum_es_ = running_sum / np.maximum(running_sum_null_mu + 2 * running_sum_null_std, 1e-9)
        running_sum_es_[-1] = float('nan')
        ax11.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            running_sum,
            c=colors[row['background']],
            linewidth=1,
        )
        ax11.plot(
            np.arange(running_sum_null_mu.shape[0]) / running_sum_null_mu.shape[0],
            running_sum_null_mu-2*running_sum_null_std,
            c=colors[row['background']],
            linestyle='dashdot',
            linewidth=1,
        )
        ax11.plot(
            np.arange(running_sum_null_mu.shape[0]) / running_sum_null_mu.shape[0],
            running_sum_null_mu,
            c=colors[row['background']],
            linestyle='dashed',
            linewidth=1,
        )
        ax11.plot(
            np.arange(running_sum_null_mu.shape[0]) / running_sum_null_mu.shape[0],
            running_sum_null_mu+2*running_sum_null_std,
            c=colors[row['background']],
            linestyle='dashdot',
            linewidth=1,
        )
        ax21.plot(
            np.arange(running_sum.shape[0]//10) / running_sum.shape[0],
            running_sum[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
        ax12.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            running_sum / row['esnull_mu'],
            c=colors[row['background']],
            linewidth=1,
        )
        ax22.plot(
            np.arange(running_sum.shape[0]//10) / running_sum.shape[0],
            (running_sum / row['esnull_mu'])[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
        ax32.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            running_sum / (running_sum_null_mu+1),
            c=colors[row['background']],
            linewidth=1,
        )
        ax42.plot(
            np.arange(running_sum.shape[0]//10) / running_sum.shape[0],
            (running_sum / (running_sum_null_mu+1))[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
        ax13.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            ((running_sum - row['esnull_mu'])/row['esnull_std']),
            c=colors[row['background']],
            linewidth=1,
        )
        ax23.plot(
            np.arange(running_sum.shape[0]//10) / running_sum.shape[0],
            ((running_sum - row['esnull_mu'])/row['esnull_std'])[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
        ax14.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            running_sum / (row['esnull_mu']+2*row['esnull_std']),
            c=colors[row['background']],
            linewidth=1,
        )
        ax24.plot(
            np.arange(running_sum.shape[0]//10) / running_sum.shape[0],
            (running_sum / (row['esnull_mu']+2*row['esnull_std']))[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
        ax33.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            running_sum_es_,
            c=colors[row['background']],
            label=f"{row['background']} (NES={row['nes']:.3f}, p={row['pval']:.3f}, ZES={row['z_es']:.3f}, ZES_p={row['z_es_pval']:.3f}, KS={row['ks_stat']:.3f}, KS_p={row['ks_p']:.3f})",
            linewidth=1,
        )
        ax43.plot(
            np.arange(running_sum.shape[0]//10) / (running_sum.shape[0]),
            running_sum_es_[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
        ax34.plot(
            np.arange(running_sum.shape[0]) / running_sum.shape[0],
            running_sum_z,
            c=colors[row['background']],
            linewidth=1,
        )
        ax44.plot(
            np.arange(running_sum.shape[0]//10) / (running_sum.shape[0]),
            running_sum_z[:running_sum.shape[0]//10],
            c=colors[row['background']],
            linewidth=1,
        )
    ax11.set_title('ES')
    ax12.set_title(r'$NES = \frac{ES}{\mu_{ESnull}}$')
    ax13.set_title(r'$NES^\prime = \frac{ES}{\mu_{ESnull} + 2 * \sigma_{ESnull}}$')
    ax14.set_title(r'$ZES = \frac{ES_x - \mu_{ESnull}}{\sigma_{ESnull}}$')

    ax32.set_title(r'$NES_x = \frac{ES_x}{\mu_{ESnull,x}+1}$')
    ax33.set_title(r'$NES^\prime_x = \frac{ES_x}{\mu_{ESnull,x} + 2 * \sigma_{ESnull,x}}$')
    ax34.set_title(r'$ZES_x = \frac{ES_x - \mu_{ESnull,x}}{\sigma_{ESnull,x}}$')
    # ax2.set_title(r'$zNES = \frac{ES - \mu_{ESnull}}{\sigma_{ESnull}}$')
    title = f"{label}"
    l = fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.05), ncol=1)
    fig.suptitle(title)
    plt.savefig(f"data/2023-04-28-{re.sub(r'[^A-Za-z]', '', title.lower())}.pdf", bbox_extra_artists=(l,), bbox_inches='tight')
    plt.show()


# %%
bgs = [
    'morgan4f',
    'physio_props',
    'L1000 Similarity',
    'mean_rank[morgan4f,physio_props,L1000 Similarity]',
    'top_rank[morgan4f,physio_props,L1000 Similarity]',
    'linear_combined[morgan4f,physio_props,L1000 Similarity]',
]
colors = dict(zip(bgs, glasbey.create_palette(len(bgs), lightness_bounds=(0, 64))))

tasks = {}
for label, records in all_results.groupby('label'):
    if label not in {'FDA Category D & X', 'Placenta Crossing'}: continue
    for _, row in records.iterrows():
        if row['background'] not in colors: continue
        if label not in tasks: tasks[label] = []
        tasks[label].append((row['background'], row))

for label, ts in tqdm(tasks.items()):
    fig, ax = plt.subplots(2, 1, figsize=(6, 6), height_ratios=[3, 1])
    labels = []
    ax = np.array(ax)
    for i, (background, r) in enumerate(tqdm(ts)):
        y = pd.Series(r['running_sum']).diff()>0
        y = y[y].index.values / r['running_sum'].shape[0]
        w = np.ones(r['running_sum'].shape[0])
        ax[0].plot(
            np.linspace(0, 1, r['running_sum'].shape[0]),
            r['running_sum'] / r['esnull_mu'],
            c=colors[background],
            label=f"{background} (NES={r['nes']:0.2f})",
        )
        ax[1].vlines(
            x=y,
            ymin=i,
            ymax=i+1,
            colors=colors[background],
        )
        labels.append(background)
        # ax[2].plot(np.arange(w.shape[0]), w, c=colors[background])

    fig.suptitle(label)
    ax[0].set_ylabel('NES')
    # ax[0].set_ylabel(r'$NES=ES/ES∅_\mu$')
    ax[0].set_xticks([])
    ax[0].set_xlim(0, 1)
    ax[1].set_yticks(np.arange(len(labels))+0.5)
    ax[1].set_yticklabels(labels)
    # ax[1].set_xticks([])
    ax[1].set_xlim(0, 1)
    # ax[2].set_xlim(0, 1)
    l = ax[0].legend(ncol=1, loc='center right', bbox_to_anchor=(-0.1, 0.5))
    # l = fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.05), ncol=1, fontsize=16)
    plt.savefig(f"data/fig-4-{re.sub(r'[^A-Za-z]', '', label.lower())}.pdf", bbox_extra_artists=(l,), bbox_inches='tight')
    plt.show()


