import api.sample_set
import os, glob, yaml

class Dataset:

    df = '' # name on file system
    fs_prefix = '' # prefix on file_system
    sample_sets = {} # Dict of sample sets, one for each preprocessing

    def __init__(self, df):
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

        for p in preprocs:
            samples = api.sample_set.SampleSet()
            samples.add_samples(self.fs_prefix, self.df, p)
            samples = samples.samples_pd[['fs_name', 'reads']].to_dict(orient='records')
            preprocessing.update({p:samples})

        self.sample_sets = preprocessing

    def __str__(self):
        return self.df + '\n' + self.fs_prefix +'\n' + str(self.sample_sets)

    def to_dict(self):
        return {
            'df': self.df,
            'fs_prefix': self.fs_prefix,
            'preprocs': self.sample_sets
        }
