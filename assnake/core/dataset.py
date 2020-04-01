import os, glob, yaml, time
import pandas as pd
from assnake.api.loaders import load_df_from_db, load_sample, load_sample_set

from assnake.utils import load_config_file,load_wc_config
from assnake.viz import plot_reads_count_change
import click

class Dataset:

    df = '' # name on file system
    fs_prefix = '' # prefix on file_system
    full_path = ''
    sample_sets = {} # Dict of sample sets, one for each preprocessing

    sources = None
    biospecimens = None
    mg_samples = None


    def __init__(self, df):
        config = load_config_file()
        wc_config = load_wc_config()
        info = load_df_from_db(df, include_preprocs = True)

        self.df =  info['df']
        self.fs_prefix =  info['fs_prefix']
        self.full_path = os.path.join(self.fs_prefix, self.df)

        preprocs = info['preprocs']
        preprocessing = {}
        for p in preprocs:
            samples = load_sample_set(wc_config, self.fs_prefix, self.df, p)
            if len(samples) > 0:
                samples = samples[['preproc', 'df', 'fs_prefix', 'df_sample', 'reads']]
                preprocessing.update({p:samples})
            

        self.sample_sets = preprocessing
        self.sample_containers = pd.concat(self.sample_sets.values())
        self.self_reads_info = self.sample_containers.pivot(index='df_sample', columns='preproc', values='reads')
  

    def plot_reads_loss(self, preprocs = [], sort = 'raw'):
        if len(preprocs) == 0: 
            preprocs = list(self.self_reads_info.columns)
        plot_reads_count_change(self.self_reads_info[preprocs].copy(), preprocs = preprocs, sort = sort, plot=True)

    def __str__(self):
        return self.df + '\n' + self.fs_prefix +'\n' + str(self.sample_sets)

    def __repr__(self):
        preprocessing_info = ''
        preprocs = list(self.sample_sets.keys())
        for preproc in preprocs:
            preprocessing_info = preprocessing_info + 'Samples in ' + preproc + ' - ' + str(len(self.sample_sets[preproc])) + '\n'
        return 'Dataset name: ' + self.df + '\n' + \
            'Filesystem prefix: ' + self.fs_prefix +'\n' + \
            'Full path: ' + os.path.join(self.fs_prefix, self.df) + '\n' + preprocessing_info

    def to_dict(self):
        preprocs = {}
        for ss in self.sample_sets:
            preprocs.update({ss : self.sample_sets[ss].to_dict(orient='records')})
        return {
            'df': self.df,
            'fs_prefix': self.fs_prefix,
            'preprocs': preprocs
        }
