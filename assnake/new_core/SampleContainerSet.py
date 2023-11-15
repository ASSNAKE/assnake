import os, glob
from datetime import datetime

import pandas as pd
from assnake.new_core.Dataset import Dataset

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
    
    def to_tsv(self, file_path, overwrite=False):
        """
        Exports the data in the SampleContainerSet to a TSV file.

        Args:
            file_path (str): The path to the file where the data should be written.
            overwrite (bool, optional): Whether to overwrite the file if it already exists.
        """
        if os.path.exists(file_path) and not overwrite:
            raise FileExistsError(f"The file '{file_path}' already exists. Use overwrite=True to overwrite it.")

        # Create a DataFrame from the SampleContainer objects
        data = []
        for container in self.sample_containers:
            container_data = {
                'df_sample': container.sample_id,
                'df': container.dataset_name,
                'fs_prefix': container.fs_prefix,
                'preproc': container.preprocessing,
            }
            data.append(container_data)

        df = pd.DataFrame(data)

        df.to_csv(file_path, sep='\t', index=False)

    def __str__(self):
        return f"SampleContainerSet for {self.dataset.dataset_name} with preprocessing {self.preproc}, {len(self.sample_containers)} containers"

    def __repr__(self):
        return f"<SampleContainerSet: {len(self.sample_containers)} containers>"
    