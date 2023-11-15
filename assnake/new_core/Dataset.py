import os, glob, datetime, yaml

import pandas as pd

from assnake.core.config import load_wc_config, read_assnake_instance_config
from assnake.new_core.Sample import Sample

class Dataset:
    """
    Represents a dataset containing multiple biological samples.
    """

    def __init__(self, dataset_name: str):
        """
        Initializes the Dataset with a name and the path to its directory.
        """
        # Load configuration and dataset information
        instance_config = read_assnake_instance_config()
        df_info_loc = os.path.join(instance_config['assnake_db'], 'datasets', dataset_name, 'df_info.yaml')
        
        if not os.path.isfile(df_info_loc):
            raise FileNotFoundError(f'No dataset information file found for dataset: {dataset_name}')
        
        with open(df_info_loc, 'r') as stream:
            df_info = yaml.safe_load(stream)


        self.dataset_name = dataset_name
        self.dataset_directory = df_info.get('fs_prefix', '')
        self.fs_prefix = df_info.get('fs_prefix', '')
        self.full_path = os.path.join(self.dataset_directory, self.dataset_name)

        self.samples = self._load_samples()

    def _load_samples(self) -> dict:
        """
        Loads all samples present in the dataset directory.
        """
        samples = {}
        raw_dir = os.path.join(self.full_path, 'reads', 'raw')
        for read_file in os.listdir(raw_dir):
            # Assuming read files are named as 'SAMPLEID_*_R1.fastq.gz'
            if '_R1' in read_file and read_file.endswith('.fastq.gz'):
                sample_id = read_file.split('_R1')[0]
                samples[sample_id] = Sample(sample_id, self.full_path)
        return samples

    def get_sample(self, sample_id: str) -> Sample:
        """
        Retrieves a sample by its ID.
        """
        return self.samples.get(sample_id)
    



    def get_sample_set(self, preproc: str, include_samples: list = None, exclude_samples: list = []) -> dict:
        """
        Generate a dictionary that represents a set of samples based on inclusion or exclusion criteria.

        :param preproc: Preprocessing step to filter samples by.
        :param include_samples: List of sample IDs to explicitly include. If None, all samples are considered.
        :param exclude_samples: List of sample IDs to exclude from the set.
        :return: A dictionary with preprocessing as keys and lists of sample IDs as values.
        """

        return []

    def save_sample_set_to_file(self, preproc, include_samples=None, exclude_samples=None, file_name=None):
        """
        Saves the selected sample set to a file.

        :param preproc: Preprocessing step to filter samples by.
        :param include_samples: List of sample IDs to explicitly include. If None, all samples are considered.
        :param exclude_samples: List of sample IDs to exclude from the set.
        :param file_name: Custom file name for saving the sample set. Defaults to the current date.
        """
        # Create a sample set using the existing method
        sample_set = self.get_sample_set(preproc, include_samples, exclude_samples)

        # Determine the file name
        if not file_name:
            file_name = datetime.now().strftime('%Y-%m-%d_sample_set.tsv')

        # Convert the sample set to a DataFrame for easy file writing
        sample_df = pd.DataFrame.from_dict(sample_set, orient='index', columns=['Preprocessing'])
        sample_df.index.name = 'Sample'

        # Save to a TSV file
        sample_df.to_csv(file_name, sep='\t')
        print(f"Sample set saved to {file_name}")

    def __str__(self):
        return f"Dataset(name={self.dataset_name}, path={self.dataset_directory}, samples={self.samples})"
