import yaml
import os
import glob
import pandas as pd

class SampleSet:
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))

    # prefix, df, preproc, fs_name
    samples = []
    samples_df = None
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

        self.samples_df = pd.DataFrame(self.samples)

    def prepare_fastqc_list_multiqc(self, strand, set_name):
        fastqc_list = []

        for s in self.samples:
            fastqc_list.append(self.wc_config['fastqc_data_wc'].format(**s, strand=strand))

        dfs = list(set(self.samples_df['df']))
        
        if len(dfs) == 1:
            prefix = list(set(self.samples_df['prefix']))[0]
            sample_list = self.wc_config['multiqc_fatqc_wc'].format(
                df = dfs[0], 
                prefix = prefix, 
                strand = strand,
                sample_set=set_name)
            print(sample_list)
            
            multiqc_dir = os.path.dirname(sample_list)
            if not os.path.isdir(multiqc_dir):
                os.makedirs(multiqc_dir)
            with open(sample_list, 'x') as file:
                file.writelines('\n'.join(fastqc_list)) 

        return fastqc_list

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

    