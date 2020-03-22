from assnake.api.sample_set import SampleSet
import os, glob, yaml, time
import pandas as pd
from assnake.utils import load_config_file

class Dataset:

    df = '' # name on file system
    fs_prefix = '' # prefix on file_system
    full_path = ''
    sample_sets = {} # Dict of sample sets, one for each preprocessing

    sources = None
    biospecimens = None
    mg_samples = None


    def __init__(self, df):
        start = time.time()

        config = load_config_file()

        # read df info
        df_info_loc = config['assnake_db']+'/datasets/{df}/df_info.yaml'.format(df = df)
        with open(df_info_loc, 'r') as stream:
            try:
                info = yaml.load(stream, Loader=yaml.FullLoader)
                if 'df' in info:
                    self.df =  info['df']
                    self.fs_prefix =  info['fs_prefix']
                    self.full_path = os.path.join(self.fs_prefix, self.df)
            except yaml.YAMLError as exc:
                print(exc)

        reads_dir = os.path.join(self.fs_prefix, self.df, 'reads/*')
        preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
        preprocessing = {}

        end = time.time()
        # print(end - start)
        for p in preprocs:
            samples = SampleSet(self.fs_prefix, self.df, p)
            # print(samples.samples_pd)
            if len(samples.samples_pd):
                samples = samples.samples_pd[['preproc', 'df', 'fs_prefix', 'fs_name', 'reads']]
            else:
                break
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
        # print(end - start)

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
