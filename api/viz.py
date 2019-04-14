import scipy.cluster.hierarchy as shc

import plotly

plotly.offline.init_notebook_mode(connected=True)
import plotly.graph_objs as go
from matplotlib import colors as mcolors
import plotly.figure_factory as FF

from skbio.stats.composition import *
from scipy.spatial.distance import pdist, squareform

def plot_changes(tr_hm2, tr_mp2, plot=True):
    traces = []
    colors  = [
        '#1f77b4',  # muted blue
        '#ff7f0e',  # safety orange
        '#2ca02c',  # cooked asparagus green
        '#d62728',  # brick red
        '#9467bd',  # muted purple
        '#8c564b',  # chestnut brown
        '#e377c2',  # raspberry yogurt pink
        '#7f7f7f',  # middle gray
        '#bcbd22',  # curry yellow-green
        '#17becf'   # blue-teal
    ]
    for i, tr in enumerate(tr_hm2):
        color_ind = i%(len(colors))
        color = colors[color_ind]
        traces.append(go.Scatter(
            y = tr['dist'],
            x = tr_mp2[i]['x'],
            name = tr['name'],
            line=dict(dash = 'dash', color = color)
        ))
        traces.append(go.Scatter(
            y = tr_mp2[i]['dist'],
            x = tr_mp2[i]['x'],
            name = tr_mp2[i]['name'],
            line=dict(color = color)
        ))
        
    fig = go.Figure(data=traces)
    if plot:
        plotly.offline.iplot(fig)
    
    return traces

def boxplots_dist(data, meta, type_of_data='', plot = True):
    dist = pdist(data, 'euclidean')
    df_dist = pd.DataFrame(squareform(dist))
    df_dist.index=data.index
    df_dist.columns=data.index
    
    feature_name = 'source'
    df = df_dist
    ####
    features = set(meta[feature_name])
    # colors = list(mcolors.CSS4_COLORS.values())
    traces = []

    for feature in features:
        samples_for_source = list(meta.loc[meta[feature_name] == feature]['sample'])
        df_sub = df[df.index.isin(samples_for_source)]
        df_sub = df_sub[df_sub.columns.intersection(samples_for_source)]

        traces.append(go.Box(
            name = feature + ' ' + type_of_data,
            x= [feature]*len(df_sub.iloc[0][1:]),
            y=df_sub.iloc[0][1:]
        ))

    # layout = dict(title = 'MDS on multiple datasets',
    #       yaxis = dict(zeroline = False, title= 'MDS1',),
    #       xaxis = dict(zeroline = False, title= 'RMDS2',)
    #      )
    
    layout = go.Layout(
        boxmode='group'
    )
    
    fig = go.Figure(data=traces,
#                     layout=layout
                   )
    if plot:
        plotly.offline.iplot(fig)
    
    return traces

def plotly_heatmap(otu_table, title):
    dend_cols = shc.dendrogram(shc.linkage(otu_table, method='ward'), no_plot = True) 
    dend_rows = shc.dendrogram(shc.linkage(otu_table.T, method='ward'), no_plot = True) 
    
    reordered_table = otu_table.T

    cols = dend_cols['ivl']
    cols_pre = list(otu_table.T.columns)
    dictionary = {}
    for i, r in enumerate(cols_pre):
        dictionary[i]=r
    new_cols = []
    for c in cols:
        new_cols.append(dictionary[int(c)])
    reordered_table = reordered_table[new_cols]
    
    
    rows = dend_rows['ivl']
    rows_pre = list(otu_table.T.index)
    dictionary = {}
    for i, r in enumerate(rows_pre):
        dictionary.update({i: r})
    new_rows = []
    for r in rows:
        new_rows.append(dictionary[int(r)])
    reordered_table = reordered_table.reindex(new_rows)
    
    
    
    figure = FF.create_dendrogram(
        otu_table, orientation='bottom', labels=cols_pre,
        linkagefun=lambda x: shc.linkage(otu_table, 'ward', metric='euclidean')
    )
    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'
        
    dendro_side = FF.create_dendrogram(
        otu_table, orientation='right',
        labels=rows_pre,
        linkagefun=lambda x: shc.linkage(otu_table.T, 'ward', metric='euclidean')
    )
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    figure.add_traces(dendro_side['data'])
    
    heatmap = [go.Heatmap( 
            z=reordered_table.values.tolist(),
            x=figure['layout']['xaxis']['tickvals'],
            y=dendro_side['layout']['yaxis']['tickvals'],
            colorscale='Viridis')]

     
    figure.add_traces(heatmap)
    # Edit xaxis
    figure['layout']["autosize"] =  True
    figure['layout']['title'].update({'text':title}),
    figure['layout']['xaxis'].update({'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks': 'outside'
                                     })
    # Edit xaxis2
    figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       'showticklabels': False,
                                       'ticks': ''
                                       }})

    # Edit yaxis
    figure['layout']['yaxis'].update({'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': True,
                                      'ticktext':dendro_side['layout']['yaxis']['ticktext'], 
                                      'tickvals':dendro_side['layout']['yaxis']['tickvals'],
                                      'ticks': 'outside'
                                     })
    # Edit yaxis2
    figure['layout'].update({'yaxis2':{'domain':[.825, .975],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       'showticklabels': False,
                                       'ticks': ''
                                      }})

    plotly.offline.iplot(figure, config={'showLink': True})
    return figure


def viz_zeroes(otu_table):
    data_zeros = otu_table[:]==0
    data_zeros = data_zeros.astype(int)
    fig = plotly_heatmap(data_zeros.T)
    
    
def time_filled_scatter(otu_table):
    x=list(otu_table.index)

    data = []
    for col in otu_table.columns:
        data.append(dict(
            x=x,
            y=otu_table[col],
            hoverinfo='x+y',
            mode='lines',
            stackgroup='one',
            name=col.split('|')[-1]
        ))

    fig = dict(data=data)
    plotly.offline.iplot(fig, config={'showLink': True})

def plot_reads_bps(meta):
    meta = meta.sort_values('bps')
    data = [go.Bar( x=meta['fs_name'], y=meta['bps'] )]
    fig = go.Figure(data=data)
    plotly.offline.iplot(fig)
    
def plot_mds(mds, feature_name, meta):
    features = set(meta[feature_name])
    colors = list(mcolors.CSS4_COLORS.values())
    traces = []

    for feature in features:
        samples_for_source = list(meta.loc[meta[feature_name] == feature]['fs_name'])
        mds_sub = mds[mds.index.isin(samples_for_source)]
        trace = go.Scatter3d(
                    x = mds_sub[0], y = mds_sub[1], z = mds_sub[2],
                    mode = 'markers',
                    marker = dict(
                        size = 6,
                    ),
                    name = feature,
                    text = mds_sub.index)
        traces.append(trace)

    layout = dict(title = 'MDS on multiple datasets',
          yaxis = dict(zeroline = False, title= 'MDS1',),
          xaxis = dict(zeroline = False, title= 'RMDS2',)
         )

    fig = go.Figure(data=traces,
                    layout=layout)
    plotly.offline.iplot(fig)