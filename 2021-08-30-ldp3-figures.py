# %%
import h5py
import numpy as np
import pandas as pd
from umap import UMAP

# %%
def iter_chunks(X, cs=1024):
    ''' Use h5py's iter_chunks method or simulate it
    '''
    import itertools
    if hasattr(X, 'iter_chunks'):
        yield from X.iter_chunks()
    else:
        for s in itertools.product(*[
            range(d//cs+int(d%cs!=0))
            for d in X.shape
        ]):
            yield tuple(
                slice(ss*cs, min(d, (ss+1)*cs))
                for ss, d in zip(s, X.shape)
            )

def iterative_mu_std(X, axis=0, cs=1024):
    ''' Compute mean on chunk generator
    '''
    mu = np.zeros(X.shape[1-axis])
    var = np.zeros(X.shape[1-axis])
    b = tuple(slice(None) if i == 1 - axis else np.newaxis for i in range(len(X.shape)))
    N = X.shape[axis]
    for s in iter_chunks(X, cs=cs):
        x_ij = X[s]
        mu[s[1-axis]] += np.sum(x_ij, axis=axis) / N
    for s in iter_chunks(X, cs=cs):
        x_ij = X[s]
        var[s[1-axis]] += np.sum(np.square(x_ij - mu[s[1-axis]][b]), axis=axis) / N
    return np.squeeze(mu), np.squeeze(np.sqrt(var))

def test_iterative_mu_std():
    ''' Verify the chunked version gives the same results on random data
    '''
    N = np.random.normal(size=(10, 5))
    mu, std = iterative_mu_std(N, axis=0, cs=3)
    assert np.isclose(mu, N.mean(axis=0)).all()
    assert np.isclose(std, N.std(axis=0)).all()
    N = np.random.normal(size=(10, 5))
    mu, std = iterative_mu_std(N.T, axis=0)
    assert np.isclose(mu, N.mean(axis=1)).all()
    assert np.isclose(std, N.std(axis=1)).all()
    N = np.random.normal(size=(10, 5))
    mu, std = iterative_mu_std(N, axis=1, cs=5)
    assert np.isclose(mu, N.mean(axis=1)).all()
    assert np.isclose(std, N.std(axis=1)).all()

test_iterative_mu_std()

# %%
lm = pd.read_csv('data/Probes_L1000_metadata.csv', index_col=0)['pr_gene_symbol']
lm

# %%
file = 'data/cp_coeff_mat.gctx'
df = h5py.File(file, 'r')
idx = pd.Series(df['0']['META']['COL']['id'][:])
lm_mask = np.in1d(df['0']['META']['ROW']['id'].asstr(), lm.values)
genes = df['0']['META']['ROW']['id'][lm_mask]

# %%
if 'matrix_norm' not in df['0/DATA/0']:
  mu, std = iterative_mu_std(df['0/DATA/0/matrix'])
  df.close()
  df = h5py.File(file, 'r+')
  df['0/DATA/0/matrix_norm'] = df['0/DATA/0'].create_dataset_like('matrix_norm', other=df['0/DATA/0/matrix'])
  for s in df['0/DATA/0/matrix'].iter_chunks():
    df['0/DATA/0/matrix_norm'][s] = (df['0/DATA/0/matrix'][s] - mu[s[1]]) / (std[s[1]] + 0.0001)
  df.close()
  df = h5py.File(file, 'r')

# %%
# this requires A LOT of memory
umap = UMAP(
    n_components=2,
    n_neighbors=50,
    min_dist=0.05,
    low_memory=True,
)
embedded = umap.fit_transform(df['0/DATA/0/matrix_norm'])
df_embedded = pd.DataFrame(embedded, columns=['UMAP-1', 'UMAP-2'], index=idx.values)

# %%
df_embedded.to_csv('data/2021-08-30-chem-umap.tsv.gz', sep='\t', compression='gzip')
