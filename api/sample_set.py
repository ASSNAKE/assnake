import yaml
import os
import glob
import pandas as pd

class SampleSet:
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))

    # prefix, df, preproc, fs_name
    samples = []
    wc_config = None
    config = None

    wc_config_loc = os.path.join(dir_of_this_file, '../wc_config.yaml')
    with open(wc_config_loc, 'r') as stream:
        try:
            wc_config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    config_loc = os.path.join(dir_of_this_file, '../config.yml')
    with open(config_loc, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    def add_samples(self, prefix, df, preproc, samples = [], do_not_add = []):
        fs_names = [f.split('/')[-1] for f in glob.glob(self.wc_config['sample_dir_wc'].format(prefix=prefix, df=df, preproc=preproc, sample = '*'))]

        for fs_name in fs_names:
            if not fs_name in do_not_add:
                self.samples.append(dict(prefix=prefix, df=df, preproc=preproc, fs_name = fs_name))        

    def prepare_dada2_sample_list(self, set_name):
        dada2_set_dir = os.path.join(self.config['dada2_dir'], set_name)

        dada2_dicts = []
        for s in self.samples:
            dada2_dicts.append(dict(mg_sample=s['fs_name'],
            R1 = self.wc_config['fastq_gz_file_wc'].format(prefix=s['prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], strand = 'R1'), 
            R2 = self.wc_config['fastq_gz_file_wc'].format(prefix=s['prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], strand = 'R2'),
            merged = self.wc_config['dada2_merged_wc'].format(prefix=s['prefix'], df=s['df'], preproc=s['preproc'], sample = s['fs_name'], sample_set = set_name)))
        if not os.path.exists(dada2_set_dir):
            os.mkdir(dada2_set_dir)

        dada2_df = pd.DataFrame(dada2_dicts)
        dada2_df.to_csv(os.path.join(dada2_set_dir, 'samples.tsv'), sep='\t', index=False)

    