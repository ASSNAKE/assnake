import os, glob, yaml, time
import pandas as pd
from assnake.api.loaders import  load_sample, load_sample_set

from assnake.core.config import load_wc_config, read_assnake_instance_config
from assnake.viz import plot_reads_count_change
import click
from pkg_resources import iter_entry_points 

class Dataset:

    df = '' # name on file system
    fs_prefix = '' # prefix on file_system
    full_path = ''
    sample_sets = {} # Dict of sample sets, one for each preprocessing

    sources = None
    biospecimens = None
    mg_samples = None


    def __init__(self, df, include_preprocs=True):
        wc_config = load_wc_config()
        instance_config = read_assnake_instance_config()

        df_info_loc = instance_config['assnake_db']+'/datasets/{df}/df_info.yaml'.format(df = df)
        df_info = {}

        if not os.path.isfile(df_info_loc):
            raise assnake.api.loaders.InputError('NO DATASET ' + df)

        with open(df_info_loc, 'r') as stream:
            try:
                info = yaml.load(stream, Loader=yaml.FullLoader)
                if 'df' in info:
                    df_info =  info
            except yaml.YAMLError as exc:
                print(exc)

        reads_dir = os.path.join(df_info['fs_prefix'], df_info['df'], 'reads/*')
        preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
        preprocessing = {}

        self.df =  df_info['df']
        self.fs_prefix =  df_info['fs_prefix']
        self.full_path = os.path.join(self.fs_prefix, self.df)

        if include_preprocs:
            preprocessing = {}
            for p in preprocs:
                samples = load_sample_set(wc_config, self.fs_prefix, self.df, p)
                if len(samples) > 0:
                    samples = samples[['preproc', 'df', 'fs_prefix', 'df_sample', 'reads']]
                    preprocessing.update({p:samples})


        self.sample_sets = preprocessing
        
        # self.sample_containers = pd.concat(self.sample_sets.values())
        # self.self_reads_info = self.sample_containers.pivot(index='df_sample', columns='preproc', values='reads')
  
    @staticmethod
    def list_in_db():
        """
        Returns dict of dictionaries with info about datasets from fs database. Key - df name
        Mandatory fields: df, prefix
        """
        dfs = {}
        instance_config = read_assnake_instance_config()
        df_info_locs = glob.glob(instance_config['assnake_db']+'/datasets/*/df_info.yaml')
        
        for df_info in df_info_locs:
            with open(df_info, 'r') as stream:
                try:
                    info = yaml.load(stream, Loader=yaml.FullLoader)
                    if 'df' in info:
                        dfs.update({info['df']: info})
                except yaml.YAMLError as exc:
                    print(exc)
        return dfs

    def plot_reads_loss(self, preprocs = [], sort = 'raw'):
        if len(preprocs) == 0: 
            preprocs = list(self.self_reads_info.columns)
        plot_reads_count_change(self.self_reads_info[preprocs].copy(), preprocs = preprocs, sort = sort, plot=True)

    def __str__(self):
        preprocessing_info = ''
        preprocs = list(self.sample_sets.keys())
        for preproc in preprocs:
            preprocessing_info = preprocessing_info + 'Samples in ' + preproc + ' - ' + str(len(self.sample_sets[preproc])) + '\n'
        return 'Dataset name: ' + self.df + '\n' + \
            'Filesystem prefix: ' + self.fs_prefix +'\n' + \
            'Full path: ' + os.path.join(self.fs_prefix, self.df) + '\n' + preprocessing_info

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



for entry_point in iter_entry_points('assnake.plugins'):
    module_class = entry_point.load()
    for k, v in module_class.dataset_methods.items():
        setattr(Dataset, k,v)
