import scipy.cluster.hierarchy as shc

import plotly

plotly.offline.init_notebook_mode()
import plotly.graph_objs as go
from matplotlib import colors as mcolors
import plotly.figure_factory as FF
import pandas as pd
# from skbio.stats.composition import *
from scipy.spatial.distance import pdist, squareform

def plotly_heatmap(otu_table, title, cluster_rows = True, cluster_cols = True, fig_size_px = (500,500)):

    reordered_table = otu_table.copy()
    if cluster_cols:
        reordered_table = reorder_by_hclust(reordered_table)
        dendro_top = FF.create_dendrogram(
            reordered_table, orientation='bottom', 
            labels=reordered_table.columns,
            linkagefun=lambda x: shc.linkage(reordered_table.T, 'ward', metric='euclidean')
        )
        for i in range(len(dendro_top['data'])):
            dendro_top['data'][i]['yaxis'] = 'y2'
        x=dendro_top['layout']['xaxis']['tickvals']
    else:
        x = reordered_table.columns
            
    if cluster_rows:
        reordered_table = reorder_by_hclust(reordered_table, False)
        dendro_side = FF.create_dendrogram(
            reordered_table.T, orientation='right',
            labels=otu_table.index,
            linkagefun=lambda x: shc.linkage(otu_table, 'ward', metric='euclidean')
        )
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'
        y=dendro_side['layout']['yaxis']['tickvals']
    else:
        y=reordered_table.index
        
    
    heatmap = [go.Heatmap( 
            z=reordered_table.values.tolist(),
            y=y,
            x=x,
            colorscale='Viridis')]
    
    figure = go.Figure(heatmap)
    
    figure['layout']["autosize"] =  True
    figure['layout']['title'].update({'text':title})
    
    
    figure['layout'].update({'xaxis': {'domain': [0, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                       'showticklabels': True,
                                      'ticks': 'outside'
                                     }})

    figure['layout'].update({'yaxis':{'domain': [0, 1],
                                          'mirror': False,
                                          'showgrid': False,
                                          'showline': False,
                                          'zeroline': False,
                                          'showticklabels': True,
                                          'ticks': 'outside'
                                         }})
    
    
    
    if cluster_cols:
        figure.add_traces(dendro_top['data'])
        print(cluster_rows)
        figure['layout'].update({'yaxis2':{'domain':[.875, .975],
                                           'mirror': False,
                                           'showgrid': False,
                                           'showline': False,
                                           'zeroline': False,
                                           'showticklabels': False,
                                           'ticks': ''
                                          }})
        
        figure['layout'].update({'xaxis': {
                                      'ticktext':dendro_top['layout']['xaxis']['ticktext'], 
                                      'tickvals':dendro_top['layout']['xaxis']['tickvals'],
                                     }})
    if cluster_rows:
        # Edit yaxis2
        figure.add_traces(dendro_side['data'])
        figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                                   'mirror': False,
                                                   'showgrid': False,
                                                   'showline': False,
                                                   'zeroline': False,
                                                   'showticklabels': False,
                                                   'ticks': ''
                                                   }})
        
        figure['layout'].update({'yaxis':{
                                      'ticktext':dendro_side['layout']['yaxis']['ticktext'], 
                                      'tickvals':dendro_side['layout']['yaxis']['tickvals']
                                     }})
    
    
    figure['layout'].update(dict(height=fig_size_px[0], width=fig_size_px[1], showlegend=False, hovermode = 'closest'))
    plotly.offline.iplot(figure, config={'showLink': True})
    
    
    
def reorder_by_hclust(dataframe, by_columns = True):
    dataframe = dataframe.copy()
    
    if by_columns:
        cols_linkage = shc.linkage(dataframe.T, method='ward')
        cols_pre = list(dataframe.columns)
    else:
        cols_linkage = shc.linkage(dataframe, method='ward')
        cols_pre = list(dataframe.index)
        
    dend_cols = shc.dendrogram(cols_linkage, no_plot = True)
    cols = dend_cols['ivl']


    dictionary = {}
    for i, r in enumerate(cols_pre):
        dictionary[i]=r

    new_cols = []
    for c in cols:
        new_cols.append(dictionary[int(c)])
        
    if by_columns:
        dataframe = dataframe[new_cols]
    else:
        dataframe = dataframe.reindex(new_cols)
        
#     print(new_cols)

    return dataframe