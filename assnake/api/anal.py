from scipy.spatial import distance
from sklearn.manifold import MDS
from sklearn.manifold import LocallyLinearEmbedding
import pandas as pd
# from skbio.stats.composition import *
from scipy.spatial.distance import pdist, squareform

import ecopy


def locally_linear_emb(data, n_neighbors=20):
    """
    Accepts pandas df rows are samples, columns are features. Indexed by sample names.
    Returns pandas df, rows are samples, axes are columns.
    """
    model = LocallyLinearEmbedding(n_neighbors=n_neighbors, n_components=3, method='modified',
                               eigen_solver='dense')
    ll = model.fit_transform(data.values)
    ll = pd.DataFrame(ll)
    ll.index = data.index
    
    return ll

def get_time_series_dist(data, meta, type_of_data=''):
    dist = pdist(data, 'euclidean')
    df_dist = pd.DataFrame(squareform(dist))
    df_dist.index=data.index
    df_dist.columns=data.index
    
    feature_name = 'source'
    df = df_dist
    features = set(meta[feature_name])
    traces = []

    for feature in features:
        samples_for_source = list(meta.loc[meta[feature_name] == feature]['sample'])
        
        df_sub = df[df.index.isin(samples_for_source)]
        df_sub = df_sub[df_sub.columns.intersection(samples_for_source)]
        
        samples_fin_names = list(df_sub.index)
        meta_sub = meta.loc[meta['sample'].isin(samples_fin_names)]
        
        time = []
        if 'time' in meta_sub.columns:
            time = list(meta_sub['time'])[1:]
        x = list(range(0,len(df_sub.iloc[0][1:])))
        
        if len(time) > 0:
            x = time
            
        traces.append({
            'name': feature + ' ' + type_of_data,
            'dist': list(df_sub.iloc[0][1:]),
            'x': x
        })
    
    return traces


def diversity(otu_table):
    methods = ['shannon' , 'gini-simpson', 'simpson' , 'dominance',  'spRich', 'even']
    index = list(otu_table.index)
    
    methods_res = {}
    for m in methods:
        div = ecopy.diversity(otu_table, method=m, breakNA=True)
#         methods_res.append(list(div))
        methods_res.update({m:list(div)})
    
    div_pd = pd.DataFrame(methods_res, columns=methods)
    div_pd.index = otu_table.index

    return div_pd

def coda(data, rm_more_than_zeroes_percent = 0.5, viz_zeroes = False):
    data_clean = data.loc[:, (data==0).mean() < rm_more_than_zeroes_percent]
    # remove rows with all zeroes
    data_clean = data_clean.loc[(data_clean.T==0).mean() != 1, :]
    if viz_zeroes:
        # TODO plot zeroes profile
        pass
    non_zero = multiplicative_replacement(data_clean)
    non_zero = pd.DataFrame(non_zero)
    non_zero.columns = data_clean.columns
    non_zero.index = data_clean.index
    clred = pd.DataFrame(clr(non_zero))
    clred.columns = data_clean.columns
    clred.index = data_clean.index
    
    return clred
    