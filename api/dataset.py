import api.sample_set
import os, glob, yaml, time
import pandas as pd

class Dataset:

    df = '' # name on file system
    fs_prefix = '' # prefix on file_system
    sample_sets = {} # Dict of sample sets, one for each preprocessing

    sources = None
    biospecimens = None
    mg_samples = None


    def __init__(self, df):
        start = time.time()

        # read config file
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        config_loc = os.path.join(curr_dir, '../config.yml')
        with open(config_loc, 'r') as stream:
            try:
                config = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        # read df info
        df_info_loc = config['assnake_db']+'/datasets/{df}/df_info.yaml'.format(df = df)
        with open(df_info_loc, 'r') as stream:
            try:
                info = yaml.load(stream)
                if 'df' in info:
                    self.df =  info['df']
                    self.fs_prefix =  info['fs_prefix']
            except yaml.YAMLError as exc:
                print(exc)

        reads_dir = os.path.join(self.fs_prefix, self.df, 'reads/*')
        preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
        preprocessing = {}

        end = time.time()
        print(end - start)
        for p in preprocs:
            samples = api.sample_set.SampleSet()
            samples.add_samples(self.fs_prefix, self.df, p)
            samples = samples.samples_pd[['preproc', 'df', 'prefix', 'fs_name', 'reads']].to_dict(orient='records')
            preprocessing.update({p:samples})
        self.sample_sets = preprocessing
        
        mg_sample_containers = []
        for preproc, sample_set in self.sample_sets.items():
            sample_set = pd.DataFrame(sample_set)
            sample_set['preproc'] = preproc
            mg_sample_containers.append(sample_set)

        self.mg_sample_containers = pd.concat(mg_sample_containers).reset_index()   
        meta = pd.concat(mg_sample_containers).reset_index()
        # Load sources
        sources_loc = config['assnake_db']+'/datasets/{df}/sources.tsv'.format(df = df)
        if os.path.isfile(sources_loc):
            self.sources = pd.read_csv(sources_loc, sep = '\t')

        # Load biospecimens
        biospecimens_loc = config['assnake_db']+'/datasets/{df}/biospecimens.tsv'.format(df = df)
        if os.path.isfile(biospecimens_loc):
            self.biospecimens = pd.read_csv(biospecimens_loc, sep = '\t')

        # Load mg samples meta
        mg_samples_loc = config['assnake_db']+'/datasets/{df}/mg_samples.tsv'.format(df = df)
        if os.path.isfile(mg_samples_loc):
            self.mg_samples = pd.read_csv(mg_samples_loc, sep = '\t')
            meta = meta.merge(self.mg_samples, left_on = 'fs_name', right_on = 'fs_name')

        # meta = self.mg_sample_containers.merge(self.mg_samples, left_on = 'fs_name', right_on = 'fs_name')
        # meta = meta.merge(self.biospecimens, left_on='biospecimen', right_on='biospecimen')
        # meta = meta.merge(self.sources, left_on='source', right_on='source')
        self.meta = meta

        end = time.time()
        print(end - start)

    def __str__(self):
        return self.df + '\n' + self.fs_prefix +'\n' + str(self.sample_sets)

    def to_dict(self):
        return {
            'df': self.df,
            'fs_prefix': self.fs_prefix,
            'preprocs': self.sample_sets
        }
