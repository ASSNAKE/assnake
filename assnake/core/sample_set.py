import yaml, os, pkg_resources, glob
import pandas as pd
import assnake.api.loaders as loaders
import assnake.utils

class SampleSet:
    """
    Class that agglomerates samples and provides convinience functions for different tasks, 
    such as constructing list of desired results locations, or preparing lists of files for rules.
    
    Attributes:
        samples_pd (:obj:`pandas.DataFrame`): Pandas DataFrame with information about samples
    """

    # prefix, df, preproc, fs_name
    samples_pd = pd.DataFrame(columns=['df', 'fs_name', 'preproc', 'reads', 'sample'])
    reads_info = pd.DataFrame()

    # ==Do we really need this here?==
    wc_config = {}
    config = {}

    
    def __init__(self, fs_prefix, df, preproc, samples_to_add = [], do_not_add = [], pattern = ''):

        self.wc_config = assnake.utils.load_wc_config()
        self.config = assnake.utils.load_config_file()

        discovered_plugins = {
            entry_point.name: entry_point.load()
            for entry_point in pkg_resources.iter_entry_points('assnake.plugins')
        }
        for module_name, module_class in discovered_plugins.items():
            for wc_config in module_class.wc_configs:
                if wc_config is not None:
                    self.wc_config.update(wc_config)
                
        self.add_samples(fs_prefix, df, preproc, samples_to_add, do_not_add, pattern)

    def add_samples(self, fs_prefix, df, preproc, samples_to_add = [], do_not_add = [], pattern = ''):
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


    

    