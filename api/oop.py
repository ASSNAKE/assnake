import pandas as pd
import math

#from skbio.stats.composition import *

import loaders as assload
import anal as anal
import viz as viz

#import plotly
#plotly.offline.init_notebook_mode(connected=True)
#import plotly.graph_objs as go
#from scipy.spatial.distance import pdist, squareform


class Dataset:
    """
    This class holds information about Dataset (Study?).
    It also provides covinience methods for data oriented workflow.
    """

    fs_prefix = ''
    df = ''
    description = ''
    paper = ''
    doi = ''
    
    fs_samples_meta = None
    samples_meta = None
    
    mp2 = None
    hm2 = None
    
    def __init__(self, df_dict):
        self.fs_prefix = df_dict.get('fs_prefix', '')
        self.df = df_dict.get('df', '')
        self.description = df_dict.get('description', '')
        self.paper = df_dict.get('paper', '')
        self.doi = df_dict.get('doi', '')
        
        self.fs_samples_meta = assload.samples_to_pd(assload.df_full_info(self.fs_prefix, self.df))
        self.samples_meta = assload.load_samples_metadata(self.fs_prefix, self.df)
        # Combine metadata with techical data 
        self.samples_meta = self.samples_meta.merge(self.fs_samples_meta, left_on='fs_name', right_on='fs_name')
        
    def sort_samples_meta(self,sort):
        self.samples_meta=self.samples_meta.sort_values(sort)
        
    def get_time_series_dist(self, on_what = 'mp2'):
        if on_what == 'mp2':
            return anal.get_time_series_dist(self.mp2, self.samples_meta, 'mp2')
        elif on_what == 'hm2':
            return anal.get_time_series_dist(self.hm2, self.samples_meta, 'hm2')
        
    def plot_heatmaps(self):
        fig = viz.plotly_heatmap(self.hm2, 'Humann2 ' + self.df)
        fig = viz.plotly_heatmap(self.mp2, 'Metaplan2 ' + self.df)
        
        
    def load_hm2(self, db, drop = True, clr = True):
        """
        This method loads Humann2 data for samples in df. 
        """
        hm2 = assload.load_hm2(
            self.fs_prefix,
            self.samples_meta.to_dict(orient='records'), 
            dbs=db, 
            index_by='sample', 
            norm = False, 
            modifier = 'norm_unstratified')

        if drop:
            cols_to_drop = ['UNMAPPED', 'UNINTEGRATED']
            hm2 = hm2.drop(cols_to_drop, axis = 1)
            
        if clr:
            hm2 = anal.coda(hm2)
            
        self.hm2 = hm2
    
    def load_mp2(self, clr = True, rm_more_than_zeroes_percent = 0.5, level='g__'):
        mp2 = assload.load_mp2(self.fs_prefix,
                               self.samples_meta.to_dict(orient='records'), 
                               level=level, 
                               index_by='sample')
        
        if clr:
            mp2=anal.coda(mp2, rm_more_than_zeroes_percent)
        self.mp2 = mp2
        
        
    def __repr__(self):
        return str({
            'fs_prefix': self.fs_prefix,
            'df': self.df,
            'description': self.description,
            'paper': self.paper,
            'doi': self.doi
        })
    
    def plot_dist_heatmaps(self, on_data = 'mp2'):
        data = None
        if on_data == 'mp2':
            data = self.mp2
        elif on_data == 'hm2': 
            data = self.hm2

            
        traces = []
        sources = set(self.samples_meta['source'])

        for source in sources:
            samples = self.samples_meta.loc[self.samples_meta['source'] == source]
            samples = list(samples['sample'])
            data_subset = data.loc[samples]

            if len(data_subset) > 1:
                coded = anal.coda(data_subset)
                dist = pdist(coded, 'euclidean')
                df_dist = pd.DataFrame(squareform(dist))

                trace = go.Heatmap(z=df_dist.values.tolist() ,
                                   x = coded.index,
                                   y = coded.index,
                                   colorscale='Viridis')
                traces.append(trace)
        
        rows = 1
        cols = 1
        
        def_cols = 3
        
        if len(traces) > def_cols:
            rows = math.ceil(len(traces)/def_cols)
            cols = def_cols
        
        curr_row = 1
        curr_col = 1
        fig = plotly.tools.make_subplots(rows=rows, cols=def_cols)
        for i, trace in enumerate(traces):
            curr_row = math.floor(i/def_cols) + 1
            
            curr_col = (i)%def_cols
            curr_col += 1
            
            fig.append_trace(trace, curr_row, curr_col)
            fig['layout']['yaxis'+str(i+1)].update(autorange='reversed')

        plotly.offline.iplot(fig, config={'showLink': True})