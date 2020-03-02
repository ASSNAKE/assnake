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
        dir_of_this_file = os.path.dirname(os.path.abspath(__file__))

        wc_config_loc = os.path.join(dir_of_this_file, '../snake/wc_config.yaml')
        with open(wc_config_loc, 'r') as stream:
            try:
                self.wc_config = yaml.load(stream, Loader=yaml.FullLoader)
            except yaml.YAMLError as exc:
                print(exc)

        discovered_plugins = {
            entry_point.name: entry_point.load()
            for entry_point in pkg_resources.iter_entry_points('assnake.plugins')
        }
        for module_name, module_class in discovered_plugins.items():
            for wc_config in module_class.wc_configs:
                if wc_config is not None:
                    self.wc_config.update(wc_config)
        
        self.config = assnake.utils.load_config_file()
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
        self.reads_info = pd.DataFrame(self.samples_pd['reads']) # Why do we need this?

    def prepare_dada2_sample_list(self, set_name='sample_set'):
        dfs = list(set(self.samples_pd['df']))
        if len(dfs) == 1:
            fs_prefix = list(set(self.samples_pd['fs_prefix']))[0]
        print(dfs)
            
        dada2_set_dir = '{fs_prefix}/{df}/dada2/{sample_set}/'.format(fs_prefix = fs_prefix, df = dfs[0], sample_set = set_name)

        dada2_dicts = []
        for s in self.samples_pd.to_dict(orient='records'):
            dada2_dicts.append(dict(mg_sample=s['fs_name'],
            R1 = self.wc_config['fastq_gz_file_wc'].format(fs_prefix=s['fs_prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], strand = 'R1'), 
            R2 = self.wc_config['fastq_gz_file_wc'].format(fs_prefix=s['fs_prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], strand = 'R2'),
            # merged = self.wc_config['dada2_merged_wc'].format(prefix=s['fs_prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], sample_set = set_name)
            ))
        if not os.path.exists(dada2_set_dir):
            os.makedirs(dada2_set_dir, exist_ok=True)

        dada2_df = pd.DataFrame(dada2_dicts)
        if not os.path.isfile(os.path.join(dada2_set_dir, 'samples.tsv')):
            dada2_df.to_csv(os.path.join(dada2_set_dir, 'samples.tsv'), sep='\t', index=False)


    def prepare_mothur_set(self, dir_loc, set_name):
        mothur_set_dir = os.path.join(dir_loc, set_name)

        mothur_dicts = []
        for s in self.samples_pd.to_dict(orient='records'):
            mothur_dicts.append(dict(mg_sample=s['fs_name'],
                R1 = self.wc_config['fastq_file'].format(prefix=s['prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], strand = 'R1'), 
                R2 = self.wc_config['fastq_file'].format(prefix=s['prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], strand = 'R2')))
        if not os.path.exists(mothur_set_dir):
            os.makedirs(mothur_set_dir)

        mothur_df = pd.DataFrame(mothur_dicts)
        mothur_df.to_csv(os.path.join(mothur_set_dir, 'stability.files'), columns=['mg_sample', 'R1', 'R2'], sep='\t', index=False, header=False)

    def prepare_assembly_set(self, assembler, params, set_name):
        # prepare dataframe
        samples = self.samples_pd[['df', 'fs_name', 'preproc']]
        samples = samples.rename({'fs_name': 'sample'}, axis=1)

        # Here we select df and prefix where the list will be saved. What to do if thre is more than 1 df i the list?
        dfs = list(set(self.samples_pd['df']))
        prefix = list(set(self.samples_pd['prefix']))[0]

        sample_table_loc = self.wc_config['assembly_table_wc'].format(
                df = dfs[0], 
                prefix = prefix, 
                assembler = assembler,
                params = params,
                sample_set=set_name)

        sample_table_dir = os.path.dirname(sample_table_loc)
        if not os.path.isdir(sample_table_dir):
            os.makedirs(sample_table_dir)
        samples.to_csv(sample_table_loc, sep='\t', index=False)

        # for s in self.samples:
        #     fastqc_list.append(self.wc_config['fastqc_data_wc'].format(**s, strand=strand))

        # 
        
        # if len(dfs) == 1:
        #     prefix = list(set(self.samples_df['prefix']))[0]
        #     sample_list = self.wc_config['multiqc_fatqc_wc'].format(
        #         df = dfs[0], 
        #         prefix = prefix, 
        #         strand = strand,
        #         sample_set=set_name)
        #     print(sample_list)
            
        #     multiqc_dir = os.path.dirname(sample_list)
        #     if not os.path.isdir(multiqc_dir):
        #         os.makedirs(multiqc_dir)
        #     with open(sample_list, 'x') as file:
        #         file.writelines('\n'.join(fastqc_list)) 

    def get_locs_for_result(self, result, preproc='', params='def'):
        result_locs = []
        strands = ['R1', 'R2']

        if result == 'count':
            for s in self.samples_pd.to_dict(orient='records'):
                if preproc == '':
                    preprocessing = s['preproc']
                else:
                    preprocessing = preproc
                for strand in strands:
                    result_locs.append(self.wc_config['count_wc'].format(
                        fs_prefix = s['fs_prefix'].rstrip('\/'),
                        df = s['df'],
                        preproc = s['preproc'],
                        sample = s['fs_name'],
                        strand = strand
                    )) 
                    # print(result_locs)
        elif result == 'fastqc':
            for s in self.samples_pd.to_dict(orient='records'):
                if preproc == '':
                    preprocessing = s['preproc']
                else:
                    preprocessing = preproc
                for strand in strands:
                    result_locs.append(self.wc_config['fastqc_zip_wc'].format(
                        fs_prefix = s['fs_prefix'].rstrip('\/'),
                        df = s['df'],
                        preproc = preprocessing,
                        sample = s['fs_name'],
                        strand = strand
                    ))
        elif result == 'mp2':
            for s in self.samples_pd.to_dict(orient='records'):
                if preproc == '':
                    preprocessing = s['preproc']
                else:
                    preprocessing = preproc
                result_locs.append(self.wc_config['mp2_out'].format(
                    fs_prefix = s['fs_prefix'].rstrip('\/'),
                    df = s['df'],
                    preproc = preprocessing,
                    sample = s['fs_name'],
                    version = '__v2.9.12'
                ))
        elif result == 'trimmomatic':
            for s in self.samples_pd.to_dict(orient='records'):
                if preproc == '':
                    preprocessing = s['preproc']
                else:
                    preprocessing = preproc
                result_locs.append(self.wc_config['fastq_gz_file_wc'].format(
                    fs_prefix = s['fs_prefix'].rstrip('\/'),
                    df = s['df'],
                    preproc = preprocessing+'__tmtic_'+params,
                    sample = s['fs_name'],
                    strand = 'R1'
                ))
        elif result == 'dada2-filter-and-trim':
            for s in self.samples_pd.to_dict(orient='records'):
                if preproc == '':
                    preprocessing = s['preproc']
                else:
                    preprocessing = preproc
                result_locs.append(self.wc_config['dada2_fat_wc'].format(
                    prefix = s['prefix'].rstrip('\/'),
                    df = s['df'],
                    preproc = preprocessing,
                    sample = s['fs_name'],
                    params = params
                ))

        
        return result_locs

    def general_taxa(self):
        samples = self.samples_pd.to_dict(orient='records')
        taxa = pd.DataFrame(loaders.get_general_taxa_comp_krak_style(samples))
        taxa.index = taxa['sample']
        return taxa
    def __str__(self):
        print('Number of samples: ', len(self.samples_pd))

    def __repr__(self):
        return 'Number of samples: ' +  str(len(self.samples_pd))

    