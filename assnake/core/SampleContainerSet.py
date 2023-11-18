import os, glob
from datetime import datetime
import shutil
import click

from assnake.utils.fs_helpers import check_and_delete_empty_directory

import pandas as pd
from assnake.core.Dataset import Dataset
import zlib

class SampleContainerSet:
    """
    Represents a set of sample containers, typically used to manage group operations in a pipeline.
    """

    def __init__(self, dataset: Dataset, preproc: str, samples_to_add: list = None, exclude_samples: list = None):
        """
        Initializes the SampleContainerSet with a dataset, a preprocessing step, and inclusion/exclusion criteria for samples.
        """
        self.dataset = dataset
        self.preproc = preproc
        self.samples_to_add = samples_to_add if samples_to_add is not None else []
        self.exclude_samples = exclude_samples if exclude_samples is not None else []
        self.sample_containers = self._build_sample_containers()

    def _build_sample_containers(self):
        """
        Builds the internal list of sample containers based on the preprocessing step and sample inclusion/exclusion lists.
        """
        sample_containers = []
        for sample_id, sample in self.dataset.samples.items():
            if (self.samples_to_add and sample_id in self.samples_to_add) or (not self.samples_to_add and sample_id not in self.exclude_samples):
                container = sample.get_container_by_preprocessing(self.preproc)
                if container:
                    sample_containers.append(container)
        return sample_containers

    @classmethod
    def create_from_kwargs(cls, df: str, preproc: str, samples_to_add: list = None, exclude_samples: list = None):
        """
        Class method to create a SampleContainerSet instance from a dataset name and preprocessing details.
        """
        dataset = Dataset(df)  # Assume a Dataset constructor that can initialize from name
        return cls(dataset, preproc, samples_to_add, exclude_samples)
    
    def rename_sample_set_with_hash(self, sample_set_file_path):
        """
        Renames the sample_set directory and the TSV file to include a hash suffix.

        Args:
            sample_set_file_path (str): Path to the sample_set.tsv file.
        """

        def compute_hash(file_path):
            """Compute crc32 hash of the TSV file content."""
            with open(file_path, 'rb') as f:
                file_content = f.read()
                hash_crc32 = zlib.crc32(file_content)
                return format(hash_crc32, 'x')

        # Extract directory and file components
        sample_set_dir = os.path.dirname(sample_set_file_path)
        sample_set_file_name = os.path.basename(sample_set_file_path)

        # Compute hash for the TSV file
        hash_suffix = compute_hash(sample_set_file_path)

        # Construct new directory and file names with hash
        new_dir_name = f"{os.path.basename(sample_set_dir)}_{hash_suffix}"
        new_dir_path = os.path.join(os.path.dirname(sample_set_dir), new_dir_name)

        # Update file name with hash suffix
        new_file_name = f"{os.path.splitext(sample_set_file_name)[0]}_{hash_suffix}.tsv"
        new_file_path = os.path.join(new_dir_path, new_file_name)

        # Create new directory if it doesn't exist
        if not os.path.exists(new_dir_path):
            os.makedirs(new_dir_path, exist_ok=True)

        # Move and rename the TSV file
        shutil.move(sample_set_file_path, new_file_path)

        check_and_delete_empty_directory(sample_set_dir)

        return new_dir_path, new_file_path

    def to_tsv(self, file_path, overwrite=False):
        """
        Exports the data in the SampleContainerSet to a TSV file and renames the dir to include a hash suffix.

        Args:
            file_path (str): The path to the file where the data should be written.
            overwrite (bool, optional): Whether to overwrite the file if it already exists.

        Returns:
            str: The path of the new file with the hash suffix.
        """
        if os.path.exists(file_path) and not overwrite:
            raise FileExistsError(f"The file '{file_path}' already exists. Use overwrite=True to overwrite it.")
        elif os.path.exists(file_path) and overwrite:
            click.secho('Overwritten', fg='yellow')

        self.tsv_path = file_path
        os.makedirs(os.path.dirname(self.tsv_path), exist_ok=True)
            
        # Create a DataFrame from the SampleContainer objects
        data = [{'df_sample': container.sample_id,
                 'df': container.dataset_name,
                 'fs_prefix': container.fs_prefix,
                 'preproc': container.preprocessing}
                for container in self.sample_containers]

        df = pd.DataFrame(data)
        df.to_csv(file_path, sep='\t', index=False)

        new_dir_path, new_file_path = self.rename_sample_set_with_hash(file_path)

        # Return the new file path with the hash
        return new_file_path

    def __str__(self):
        return f"SampleContainerSet for {self.dataset.dataset_name} with preprocessing {self.preproc}, {len(self.sample_containers)} containers"

    def __repr__(self):
        return f"<SampleContainerSet: {len(self.sample_containers)} containers>"
    