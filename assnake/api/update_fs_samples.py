from tempfile import NamedTemporaryFile
import shutil
import csv
import os
from assnake import Dataset
from assnake import utils
from assnake.api import fs_helpers, loaders
import pandas as pd
fields = ['sample', 'sequencing_run']


def update_fs_samples_csv(dataset):
    '''
    Scans dataset folder and updates fs_samples.tsv if necessary
    
    :param dataset: Name of the dataset
    :return: Returns sample dict in loc
    
    '''
    fs_samples_tsv_loc = '{config}/datasets/{df}/fs_samples.tsv'.format(config=utils.load_config_file()['assnake_db'], df=dataset)
    df = Dataset(dataset)
    fs_samples_pd = pd.concat(list(df.sample_sets.values()))
    fs_samples_pd['Final'] = 'never_set'
    fs_samples_pd.to_csv(fs_samples_tsv_loc, sep='\t', index = False)

    return True