from assnake.core.sample_set import SampleSet
import os, glob, yaml, time
import pandas as pd
from assnake.api.loaders import load_df_from_db
from assnake.utils import load_config_file
from assnake.viz import plot_reads_count_change


def load_sample_set(df, preproc, fs_prefix = '', samples_to_add = [], exclude_samples = [], pattern = ''):
    '''
    This function is used to add samples into the SampleSet.

    Args:
        fs_prefix: Prefix of the dataset on filesystem
        df: Name of the dataset
        preproc: Preprocessing you want to use
        samples_to_add: List of sample names to add
        do_not_add: list of sample names NOT to add
        pattern: sample names must match this glob pattern to be included. 
    '''
    df_loaded = assnake.api.loaders.load_df_from_db(df)

    # Now for the meta column stuff
    meta_loc = os.path.join(df_loaded['fs_prefix'], df_loaded['df'], 'mg_samples.tsv')
    if os.path.isfile(meta_loc):
        meta = pd.read_csv(meta_loc, sep = '\t')
        if meta_column is not None:
            if column_value is not None:
                subset_by_col_value = meta.loc[meta[meta_column] == column_value]
                if len(subset_by_col_value) > 0:
                    samples_to_add = list(subset_by_col_value['sample_name'])
            else:
                print('dfg')


    sample_set = assnake.SampleSet(df_loaded['fs_prefix'], df_loaded['df'], preproc, samples_to_add=samples_to_add)
    if len(exclude_samples) > 0 :  
        sample_set.samples_pd = sample_set.samples_pd.loc[~sample_set.samples_pd['fs_name'].isin(exclude_samples), ]

    click.echo(tabulate(sample_set.samples_pd[['fs_name', 'reads', 'preproc']].sort_values('reads'), headers='keys', tablefmt='fancy_grid'))

    # construct sample set name for fs
    if meta_column is None and column_value is None:
        curr_date = datetime.datetime.now()
        def_name = '{month}{year}'.format(month=curr_date.strftime("%b"), year=curr_date.strftime("%y"))
        sample_set_name = def_name
    else:
        sample_set_name = meta_column + '__' + column_value

    samples = []
    fastq_gz_file_loc = self.wc_config['fastq_gz_file_wc'].format(
        fs_prefix=fs_prefix, df=df, preproc=preproc, 
        strand='R1',sample = '*')
    fs_names = [f.split('/')[-1].split('.')[0].replace('_R1', '') for f in glob.glob(fastq_gz_file_loc)]

    if pattern != '':
        fs_names = [f.split('/')[-1] for f in 
        glob.glob(self.wc_config['sample_dir_wc'].format(fs_prefix=fs_prefix, df=df, preproc=preproc, sample = pattern))]

    sample_dir_wc = self.wc_config['sample_dir_wc']
    fastq_gz_file_wc = self.wc_config['fastq_gz_file_wc']
    count_wc = self.wc_config['count_wc']

    fs_names = list(set(fs_names) - set(do_not_add))
    if len(samples_to_add) > 0:
        fs_names = fs_names and samples_to_add

    samples = [loaders.load_sample(fs_prefix, df, preproc, fs_name,
                    sample_dir_wc = sample_dir_wc, fastq_gz_file_wc = fastq_gz_file_wc, 
                    count_wc=count_wc) for fs_name in fs_names]
    
    samples_pd = pd.DataFrame(samples)
    # samples_pd.index = samples_pd['fs_name'] + ':' + samples_pd['preproc'] # We can debate on this

    self.samples_pd = pd.concat([self.samples_pd, samples_pd], sort=True)



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
        info = load_df_from_db(df, include_preprocs = True)

        self.df =  info['df']
        self.fs_prefix =  info['fs_prefix']
        self.full_path = os.path.join(self.fs_prefix, self.df)

        preprocs = info['preprocs']
        preprocessing = {}

        for p in preprocs:
            samples = SampleSet(self.fs_prefix, self.df, p)
            if len(samples.samples_pd):
                samples = samples.samples_pd[['preproc', 'df', 'fs_prefix', 'fs_name', 'reads']]
            else:
                break
            preprocessing.update({p:samples})

        self.sample_sets = preprocessing
        self.sample_containers = pd.concat(self.sample_sets.values())
        self.self_reads_info = self.sample_containers.pivot(index='fs_name', columns='preproc', values='reads')
  

    def plot_reads_loss(self, preprocs = [], sort = 'raw'):
        # preprocs = list(self.self_reads_info.columns)
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
