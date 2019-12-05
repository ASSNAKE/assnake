import pandas as pd
import math, os

#from skbio.stats.composition import *

import loaders as assload
import anal as anal
import viz as viz

import glob

import yaml


#import plotly
#plotly.offline.init_notebook_mode(connected=True)
#import plotly.graph_objs as go
#from scipy.spatial.distance import pdist, squareform

config = None
dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
cofig_loc = os.path.join(dir_of_this_file, '../config.yml')
with open(cofig_loc, 'r') as stream:
    try:
        config = yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)

ASSNAKE_DB = config.get('assnake_db')

SOURCES_META_LOC = '{assnake_db}/datasets/{df}/sources.tsv'
BIOSPECIMENS_META_LOC = '{assnake_db}/datasets/{df}/biospecimens.tsv'
MG_SAMPLES_META_LOC = '{assnake_db}/datasets/{df}/mg_samples.tsv'

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
    
    def __init__(self, df_dict, preproc = 'longest'):
        self.fs_prefix = df_dict.get('fs_prefix', '')
        self.df = df_dict.get('df', '')
        self.description = df_dict.get('description', '')
        self.paper = df_dict.get('paper', '')
        self.doi = df_dict.get('doi', '')

        self.fs_samples_meta = assload.samples_to_pd(assload.df_full_info(self.fs_prefix, self.df, preproc))


        sources_loc = SOURCES_META_LOC.format(assnake_db=ASSNAKE_DB, df=self.df)
        meta_loc = BIOSPECIMENS_META_LOC.format(assnake_db=ASSNAKE_DB, df=self.df)
        mg_samples_loc = MG_SAMPLES_META_LOC.format(assnake_db=ASSNAKE_DB, df=self.df)

        meta = pd.read_csv(meta_loc, sep='\t')
        if os.path.isfile(sources_loc):
            sources = pd.read_csv(sources_loc, sep='\t')
            meta = meta.merge(sources, left_on='source', right_on='source', suffixes=('_biospecimen', '_source'))

        mg_samples_meta = pd.read_csv(mg_samples_loc, sep='\t')
        meta = meta.merge(mg_samples_meta, left_on='biospecimen', right_on='biospecimen' )

        self.samples_meta = meta

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
            index_by='fs_name', 
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
                               index_by='fs_name')
        
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


class MagCollection:
    '''
    Wrapper class for working with collections of MAGs. 
    '''
    dfs = ''
    preprocs = ''
    samples = ''

    bins = []

    bins_wc = '/data5/bio/databases/fna/assembly/mh__def/{dfs}/{samples}/imp__tmtic_def1/{collection}/bin_by_bin/{binn}'
    taxa_wc = '/data5/bio/databases/fna/assembly/mh__def/{dfs}/{samples}/imp__tmtic_def1/{collection}/bin_by_bin/{binn}/{binn}-bin_taxonomy.tab'
    summary_wc = '/data5/bio/databases/fna/assembly/mh__def/{dfs}/{samples}/imp__tmtic_def1/{collection}/bins_summary.txt'
    checkm_wc = '/data5/bio/databases/fna/assembly/{assembler}/{dfs}/{samples}/imp__tmtic_def1/{collection}/checkm/storage/bin_stats_ext.tsv'

    def __init__(self, dfs, preprocs, samples, collection, assembler):
        self.dfs = dfs
        self.preprocs = preprocs,
        self.samples = samples
        self.collection = collection
        self.assembler = assembler
        bins  = [r.split('/')[-1] for r 
             in glob.glob(self.bins_wc.format(binn = '*', dfs = dfs, samples = self.samples, collection=collection))]

        mags = []
        for b in bins:
            try:
                taxa = pd.read_csv(self.taxa_wc.format(binn = b, dfs = dfs, samples = self.samples, collection=collection), header=None, sep='\t')
                taxa = taxa.fillna('Unknown')
                dd = {"Bin": b,
                'Taxa': list(taxa[0])[0].split('-')[0] + '__' + list(taxa[1])[0]
                }
                mags.append(dd)
            except:
                pass
                print("Can't load: ", self.taxa_wc.format(binn = b, dfs = dfs, samples = self.samples, collection=collection))
            
            # mags.append(dd)
        self.bins = pd.DataFrame(mags)

        
        try:
            self.summary = pd.read_csv(self.summary_wc.format(samples = self.samples, dfs = dfs, collection=collection), sep='\t')
            self.checkm = self.get_bins()
            subcheckm = self.checkm[['Bin', 'Completeness', 'Contamination', 'marker lineage']]
            self.summary = self.summary.merge(subcheckm, left_on = 'bins', right_on = 'Bin')
            self.summary = self.summary.drop(['Bin'], axis=1)
            self.summary = self.summary.merge(self.bins, left_on = 'bins', right_on = 'Bin')
            self.summary = self.summary.drop(['bins'], axis=1)
        except:
            pass
        

    def filter_by_comp_cont(self, completeness, contamination):
        summ = self.summary.loc[self.summary['Completeness'] > completeness]
        summ = summ.loc[summ['Contamination'] < contamination]
        return summ

    def __repr__(self):
        return str({
            'dfs': self.dfs,
            'samples': self.samples,
            'preprocs': self.preprocs,
            'bins': self.bins,
        })

    def get_bins(self):
        bins = []
        with open(self.checkm_wc.format(dfs = self.dfs, samples = self.samples, collection=self.collection,assembler=self.assembler), 'r') as checkm: 
            for line in checkm:
                dic=eval(line.strip().split("\t")[1])
                dic.update({'Bin': line.strip().split("\t")[0]})
                bins.append(dic)
        bbb = pd.DataFrame(bins)
        return bbb


class Mag:
    def __init__(self, dfs, preprocs, samples, collection, binn):
        names_wc = '/data5/bio/databases/fna/assembly/mh__def/FHM/{samples}/imp__tmtic_def1/{collection}/bin_by_bin/{binn}/{binn}-contigs.names'

        with open(names_wc.format(binn = binn, samples = samples, collection=collection), 'r') as names:
            self.contigs = [c.strip() for c in names.readlines()]

