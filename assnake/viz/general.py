import scipy.cluster.hierarchy as shc

import plotly

# plotly.offline.init_notebook_mode()
import plotly.graph_objs as go
from matplotlib import colors as mcolors
import plotly.figure_factory as FF
import pandas as pd
# from skbio.stats.composition import *
from scipy.spatial.distance import pdist, squareform
    
def plot_centr(centr, subtitle=''):
    trace4 = go.Bar(
        y=list(centr.index),
        x=list(centr['uncl']/centr['total']),
        name='uncl',
        orientation = 'h',
    )
    trace2 = go.Bar(
        y=list(centr.index),
        x=list(centr['bacteria']/centr['total']),
        name='bacteria',
        orientation = 'h',
    )
    trace3 = go.Bar(
        y=list(centr.index),
        x=list(centr['other']/centr['total']),
        name='other',
        orientation = 'h',
    )
    trace1 = go.Bar(
        y=list(centr.index),
        x=list(centr['homo']/centr['total']),
        name='HUMAN',
        orientation = 'h',
    )
    trace5 = go.Bar(
        y=list(centr.index),
        x=list(centr['vir']/centr['total']),
        name='vir',
        orientation = 'h'
    )
    trace6 = go.Bar(
        y=list(centr.index),
        x=list(centr['archaea']/centr['total']),
        name='archaea',
        orientation = 'h'
    )

    data = [trace4, trace2, trace3, trace1, trace5, trace6]
    layout = go.Layout(
        barmode='stack',
        title = 'General taxa composition; ' + subtitle 
    )

    fig = go.Figure(data=data, layout=layout)
    if False: plotly.offline.iplot(fig)

def plot_mds(mds, feature_name, meta, title='MDS', select_by='fs_name'):
    features = set(meta[feature_name])
    colors = list(mcolors.CSS4_COLORS.values())
    traces = []

    for feature in features:
        samples_for_source = list(meta.loc[meta[feature_name] == feature][select_by])
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

    layout = dict(title = title,
          yaxis = dict(zeroline = False, title= 'MDS1',),
          xaxis = dict(zeroline = False, title= 'RMDS2',)
         )

    fig = go.Figure(data=traces,
                    layout=layout)
    if False: plotly.offline.iplot(fig)



def plot_reads_count_change(read_table, preprocs, sort, title = 'Reads number', trace_names = [], number_index_hack=False, plot=False):
    """
    This function plots bar chart with reads lost on each preprocessing step.

    Args:
        read_table (:obj:`pandas.DataFrame`): DataFrame with columns with number of reads for each sample on every preprocessing step
        preprocs (list[str]): List of preprocessings in called order
        sort (list[str]): Sort by this variable
    """
    read_table.index = read_table.index.map(str) +'_'
    read_table = read_table.sort_values(sort, ascending=False)

    if len(trace_names) == 0:
        trace_names = preprocs
        
    prev_preproc = preprocs[0]
    traces = []
    for p in preprocs[1:]:
        read_table['diff_'+prev_preproc+'--'+p] = read_table[prev_preproc] - read_table[p]
        traces.append(go.Bar(x=read_table.index, y=read_table['diff_'+prev_preproc+'--'+p], name='diff_'+prev_preproc+'--'+p))
        prev_preproc = p

    
    traces.append(go.Bar( x=read_table.index, y=read_table[preprocs[-1]], name=preprocs[-1]))

    traces.reverse()
    
    layout = go.Layout( barmode='stack', margin=go.layout.Margin( b=100 ), width=1800, title = title)
    fig = go.Figure(data=traces, layout=layout)
    if plot: plotly.offline.iplot(fig)