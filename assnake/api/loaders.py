import os
import glob

import pandas as pd
import numpy as np
import yaml
import assnake
from assnake.core.config import read_assnake_instance_config, load_wc_config
from assnake.utils.general import bytes2human


# TODO this goes to Exceptions 
class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

# TODO this goes to assnake-core-preprocessing
def load_count(fs_prefix, df, preproc, df_sample, report_bps=False, verbose=False, count_wc=''):
    """
    Loads information about read and bp count in paired-end sample.
    """
    strands = ['R1', 'R2']
    
    count_loc1 = count_wc.format(fs_prefix=fs_prefix, df=df, preproc=preproc, df_sample=df_sample, strand=strands[0])
    count_loc2 = count_wc.format(fs_prefix=fs_prefix, df=df, preproc=preproc, df_sample=df_sample, strand=strands[1])
    
    count_dict = {'reads': -1}

    try:
        with open(count_loc1, 'r') as f:
            line = f.readline().rstrip()
            count_dict['reads'] += int(line.split(' ')[0]) + 1
            if report_bps:
                count_dict.update({'bps':int(line.split(' ')[1])})
        # with open(count_loc2, 'r') as f:
        #     line = f.readline().rstrip()
        #     count_dict['reads'] += int(line.split(' ')[0])
        #     if report_bps:
        #         count_dict['bps'] += int(line.split(' ')[1])
    except:
        if verbose: 
            print('error loading counts: ', sample)
        if report_bps:
            count_dict.update({'bps': -1})
        return count_dict
    
    return count_dict
        

def load_sample(fs_prefix, df, preproc, df_sample,
                report_bps=False, report_size=False, verbose=False,
                sample_dir_wc = '', fastq_gz_file_wc = '', count_wc=''):
    '''
    Loads all necessary info about given sample from file system.
    '''
    # Init start values
    sample_dict = {}
    strands = ['R1', 'R2']
    
    final_preproc = ''
    size = 0
    containers = []
    preprocs = []


    # Now select what preprocessing we want to use
    if preproc == 'longest':
        preprocs = [p.split('/')[-2] for p in 
                    glob.glob(fastq_gz_file_wc.format(fs_prefix=fs_prefix, df=df, preproc='*', df_sample=df_sample, strand='R1'))]
    else:
        preprocs = [p.split('/')[-2] for p in 
                    glob.glob(fastq_gz_file_wc.format(fs_prefix=fs_prefix, df=df, preproc=preproc, df_sample=df_sample, strand='R1'))]
    for p in preprocs:
        r1 = fastq_gz_file_wc.format(fs_prefix=fs_prefix, df=df, preproc=p, df_sample=df_sample, strand=strands[0])
        r2 = fastq_gz_file_wc.format(fs_prefix=fs_prefix, df=df, preproc=p, df_sample=df_sample, strand=strands[1])

        if os.path.isfile(r1) and os.path.isfile(r2):
            containers.append(p)
            if len(p) > len(final_preproc):
                final_preproc = p
                if report_size:
                    size = os.path.getsize(r1) + os.path.getsize(r2)
                    sample_dict.update({'size': bytes2human(size, symbols='iec'), 'bytes': size})
        elif os.path.isfile(r1):
            containers.append(p)
            if len(p) > len(final_preproc):
                final_preproc = p
                if report_size:
                    size = os.path.getsize(r1) 
                    sample_dict.update({'size': bytes2human(size, symbols='iec'), 'bytes': size})
        else:
            click.secho('There are no reads in a path: %s'%r1, fg='red')
            exit()
    return {'df':df, 
            'df_sample':df_sample, 
            'df_sample':df_sample,  
            'preproc':final_preproc, 
            'fs_prefix': fs_prefix,
            #'preprocs':containers, 
            **load_count(fs_prefix, df, final_preproc, df_sample, verbose, count_wc=count_wc)}


def load_sample_set(wc_config, fs_prefix, df, preproc, samples_to_add = [], do_not_add = [], pattern = '*'):
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
    if wc_config is None:
        wc_config = load_wc_config()

    samples = []
    fastq_gz_file_loc = wc_config['fastq_gz_file_wc'].format(
        fs_prefix=fs_prefix, df=df, preproc=preproc, 
        strand='R1', df_sample = pattern)
    
    df_samples = [f.split('/')[-1].split('.')[0].replace('_R1', '') for f in glob.glob(fastq_gz_file_loc)]
    sample_dir_wc    = wc_config['sample_dir_wc']
    fastq_gz_file_wc = wc_config['fastq_gz_file_wc']
    count_wc         = wc_config['count_wc']

    df_samples = list(set(df_samples) - set(do_not_add))

    if len(samples_to_add) > 0: df_samples = df_samples and samples_to_add

    samples = [load_sample(fs_prefix, df, preproc, df_sample,
                    sample_dir_wc = sample_dir_wc, fastq_gz_file_wc = fastq_gz_file_wc, 
                    count_wc=count_wc) for df_sample in df_samples]
    
    sample_set = pd.DataFrame(samples)
    return sample_set

fields = ['sample', 'sequencing_run']

# TODO where to put this one?
def update_fs_samples_csv(dataset):
    '''
    Scans dataset folder and updates fs_samples.tsv if necessary
    
    :param dataset: Name of the dataset
    :return: Returns sample dict in loc
    
    '''
    fs_samples_tsv_loc = '{assnake_db}/datasets/{df}/assnake_samples.tsv'.format(assnake_db=read_assnake_instance_config()['assnake_db'], df=dataset)
    df = assnake.Dataset(dataset)
    fs_samples_pd = pd.concat(list(df.sample_sets.values()))
    fs_samples_pd = df.sample_sets['raw']
    fs_samples_pd['final_preprocessing'] = 'never_set'
    fs_samples_pd.to_csv(fs_samples_tsv_loc, sep='\t', index = False)

    return True



