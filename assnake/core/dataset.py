import os
import glob
import yaml
import pandas as pd
from assnake.api.loaders import load_sample_set
from assnake.core.config import load_wc_config, read_assnake_instance_config

class Dataset:
    """
    A class to manage datasets for metagenomic analysis. It provides functionality
    to load and handle dataset preprocessing information stored in a file system.
    
    Attributes:
        df (str): Name of the dataset on the file system.
        fs_prefix (str): Filesystem prefix for the dataset path.
        full_path (str): Full file system path to the dataset.
        sample_sets (dict): A dictionary containing sample sets for each preprocessing step.
    """

    def __init__(self, df, include_preprocs=True):
        """
        Initializes a Dataset object by loading the dataset's information and preprocessing samples.

        Args:
            df (str): The name of the dataset to initialize.
            include_preprocs (bool): Flag to determine whether to load preprocessing samples.
        """
        self.df = df
        self.sample_sets = {}
        
        # Load configuration and dataset information
        instance_config = read_assnake_instance_config()
        df_info_loc = os.path.join(instance_config['assnake_db'], 'datasets', df, 'df_info.yaml')
        
        if not os.path.isfile(df_info_loc):
            raise FileNotFoundError(f'No dataset information file found for dataset: {df}')
        
        with open(df_info_loc, 'r') as stream:
            df_info = yaml.safe_load(stream)
        
        self.fs_prefix = df_info.get('fs_prefix', '')
        self.full_path = os.path.join(self.fs_prefix, self.df)

        # Load sample sets if requested
        if include_preprocs:
            self._load_sample_sets()

    def _load_sample_sets(self):
        """
        Internal method to load sample sets from the dataset's preprocessing folders.
        """
        reads_dir = os.path.join(self.full_path, 'reads', '*')
        preprocs = [os.path.basename(p) for p in glob.glob(reads_dir)]
        
        for p in preprocs:
            sample_set = load_sample_set(load_wc_config(), self.fs_prefix, self.df, p)
            if len(sample_set) > 0:
                self.sample_sets[p] = sample_set[['preproc', 'df', 'fs_prefix', 'df_sample', 'reads']]

    @staticmethod
    def list_in_db():
        """
        Static method to list all datasets available in the ASSNAKE database.

        Returns:
            dict: A dictionary of datasets with their information.
        """
        instance_config = read_assnake_instance_config()
        df_info_locs = glob.glob(instance_config['assnake_db']+'/datasets/*/df_info.yaml')
        dfs = {}
        
        for df_info_loc in df_info_locs:
            with open(df_info_loc, 'r') as stream:
                info = yaml.safe_load(stream)
                if 'df' in info:
                    dfs[info['df']] = info
        return dfs

    @staticmethod
    def delete_ds(dataset):
        """
        Static method to remove a dataset from the ASSNAKE database.

        Args:
            dataset (str): The name of the dataset to remove.

        Returns:
            tuple: A tuple containing a boolean success flag and an error message if any.
        """
        try:
            config_loc = read_assnake_instance_config()['assnake_db']
            os.remove(f"{config_loc}/datasets/{dataset}/df_info.yaml")
            return True, None
        except Exception as e:
            return False, str(e)

    def __str__(self):
        """
        String representation of the Dataset object.

        Returns:
            str: Information about the Dataset object.
        """
        preprocessing_info = '\n'.join([f'Samples in {preproc}: {len(self.sample_sets[preproc])}' 
                                        for preproc in self.sample_sets])
        return f'Dataset name: {self.df}\nPrefix: {self.fs_prefix}\nFull path: {self.full_path}\n{preprocessing_info}'

    def __repr__(self):
        """
        Technical representation of the Dataset object, similar to the string representation.
        
        Returns:
            str: Information about the Dataset object for debugging.
        """
        return self.__str__()

    def to_dict(self):
        """
        Converts the Dataset object to a dictionary format.

        Returns:
            dict: A dictionary representation of the Dataset object.
        """
        preprocs_dict = {preproc: self.sample_sets[preproc].to_dict(orient='records') 
                         for preproc in self.sample_sets}
        return {
            'df': self.df,
            'fs_prefix': self.fs_prefix,
            'preprocs': preprocs_dict
        }
