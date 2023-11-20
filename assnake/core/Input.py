from datetime import datetime
import os 
from assnake.core.SampleContainerSet import SampleContainerSet

class Input:
    def __init__(self, dataset, preprocessing=None, samples_to_add=None, exclude_samples=None, additional_input_options=None):
        self.dataset = dataset
        self.preprocessing = preprocessing
        self.samples_to_add = samples_to_add
        self.exclude_samples = exclude_samples
        self.additional_input_options = additional_input_options or {}
        if preprocessing:
            self.sample_container_set = SampleContainerSet(dataset, preprocessing, samples_to_add, exclude_samples)
        else:
            self.sample_container_set = None

    def prepare_snakemake_output_paths(self, result, preset):
        raise NotImplementedError("This method should be implemented in subclasses")


class IlluminaSampleInput(Input):
    """
    Input class for individual Illumina samples.
    """

    def __init__(self, dataset, preprocessing, samples_to_add, exclude_samples, additional_input_options):
        """
        Initializes an IlluminaSampleInput instance.

        Args:
            dataset (Dataset): The dataset associated with the input.
            preprocessing (str): The preprocessing step name.
            samples_to_add (list): A list of sample names to be included.
            exclude_samples (list): A list of sample names to be excluded.
        """
        super().__init__(dataset, preprocessing, samples_to_add, exclude_samples, additional_input_options)

    def get_input_list_for_formatting(self):
        formatting_list = []
        for container in self.sample_container_set.sample_containers:
            formatting_dict = dict(
                fs_prefix=container.fs_prefix,
                df=container.dataset_name,
                preproc=container.preprocessing,
                df_sample=container.sample_id
            )
            formatting_list.append(formatting_dict)
        return formatting_list




class IlluminaSampleSetInput(Input):
    """
    Input class for a set of Illumina samples.
    """
    def __init__(self, dataset, preprocessing, samples_to_add, exclude_samples, additional_input_options=None):
        super().__init__(dataset, preprocessing, samples_to_add, exclude_samples, additional_input_options)
        self.sample_set_name = additional_input_options.get('sample_set_name')
        print(self.additional_input_options)
        
    def create_sample_set_tsv(self):
        """
        Creates a TSV file for the sample set and returns its path.
        """
        tsv_file_path = os.path.join(
            self.dataset.full_path, 
            self.additional_input_options['subpath_for_sample_set_tsv'], 
            self.sample_set_name, 
            "sample_set.tsv")

        tsv_file_with_hash = self.sample_container_set.to_tsv(tsv_file_path, overwrite=True)
        return tsv_file_with_hash

    def get_input_list_for_formatting(self):
        self.sample_set_file = self.create_sample_set_tsv()
        sample_set_name = self.sample_set_file.split('/')[-2]

        formatting_list = [
            dict(
                fs_prefix=self.dataset.fs_prefix,
                df=self.dataset.dataset_name,
                sample_set=sample_set_name
            )
        ]
        return formatting_list
    
    def prepare_snakemake_output_paths(self, result, preset):
        """
        Prepares Snakemake output paths for the sample set.

        Args:
            result (Result): The result object associated with this input.
            preset (str): The preset name applied to the result.

        Returns:
            List[str]: A list of Snakemake target file paths.
        """
        target = result.format_target_path(sample_set=self.sample_set_name, preset=preset)
        return [target]
    

import os
import yaml
import hashlib

import re


class FeatureTableSpecificationInput(Input):
    # PRELIMINARY ROUGH IMPLEMENTATION

    """
    Input subclass for handling feature table specifications.
    """

    def __init__(self, dataset, preproc, samples_to_add, exclude_samples, additional_input_options ):
        """
        Initializes a FeatureTableSpecificationInput instance.

        Args:
            dataset (Dataset): The dataset associated with the feature table.
            sample_set (str): The name of the sample set used to create the feature table.
            feature_table_name (str): User-defined or default name for the feature table.
            source_wc_string (str): Wildcard string representing the original produced feature table.
        """
        super().__init__(dataset, None)
        self.additional_input_options = additional_input_options
        self.feature_table_name = additional_input_options['feature_table_name']
        self.source_wc_string = additional_input_options['source_wc_string']
        self.sample_set = additional_input_options['sample_set_name']
        self.metadata = self._create_metadata()
        self.feature_table_path = self._get_feature_table_path()

    def _create_metadata(self):
        """
        Creates metadata for the feature table based on the source wildcard string.
        """

        print('==========================================')
        print(self.additional_input_options)

        source_loc = self.additional_input_options.get('source_wc_string', None)
        source_loc = source_loc.format(
            df = self.dataset.dataset_name,
            fs_prefix = self.dataset.fs_prefix,
            sample_set = self.additional_input_options.get('sample_set_name'),
            learn_errors_preset = self.additional_input_options.get('learn_errors_preset'),
            core_algo_preset = self.additional_input_options.get('core_algo_preset'),
            merged_preset = self.additional_input_options.get('merged_preset'),
            nonchim_preset = self.additional_input_options.get('nonchim_preset'),
        )

        # Parse presets and other info from source_wc_string
        # parsed_metadata = self._parse_wc_string(self.source_wc_string, source_loc)

        parsed_metadata = dict(
            sample_set = self.additional_input_options.get('sample_set_name'),
            learn_errors_preset = self.additional_input_options.get('learn_errors_preset'),
            core_algo_preset = self.additional_input_options.get('core_algo_preset'),
            merged_preset = self.additional_input_options.get('merged_preset'),
            nonchim_preset = self.additional_input_options.get('nonchim_preset'),
        )
        self.metadata = parsed_metadata
        print(parsed_metadata)

        # # Additional metadata can be added here
        # parsed_metadata['sample_set'] = self.sample_set

        # Saving metadata to a YAML file
        metadata_filename = "metadata.yaml"
        metadata_filepath = os.path.join(self.dataset.full_path, 'feature_tables', self.sample_set, self.feature_table_name, metadata_filename)
        
        os.makedirs(os.path.dirname(metadata_filepath), exist_ok=True)
        with open(metadata_filepath, 'w') as metadata_file:
            yaml.dump(parsed_metadata, metadata_file)

        return parsed_metadata

    def _get_feature_table_path(self):
        """
        Computes the path for the feature table using hash of the metadata.
        """
        # Compute hash of metadata
        metadata_hash = self._compute_hash_from_metadata(self.metadata)

        # Construct feature table path
        feature_table_dir = os.path.join(self.dataset.full_path, 'feature_tables', self.sample_set, self.feature_table_name)
        feature_table_path = os.path.join(feature_table_dir, f"{metadata_hash}.rds")
        print(feature_table_path)
        return feature_table_path

    def get_input_list_for_formatting(self):

        formatting_list = [
            dict(
                fs_prefix=self.dataset.fs_prefix,
                df=self.dataset.dataset_name,
                ft_name = self.feature_table_name,  
                sample_set=self.sample_set
            )
        ]
        return formatting_list


    def _compute_hash_from_metadata(self, metadata):
        """
        Computes a hash from the metadata dictionary.
        """
        hash_obj = hashlib.md5(str(metadata).encode())
        return hash_obj.hexdigest()

    def prepare_snakemake_output_paths(self, result, preset):
        """
        Prepares Snakemake output paths for the feature table.

        Args:
            result (Result): The result object associated with this input.
            preset (str): The preset name applied to the result.

        Returns:
            List[str]: A list of Snakemake target file paths.
        """
        # Assuming that result is relevant to the feature table creation
        target = result.format_target_path(feature_table_name=self.feature_table_name, sample_set=self.sample_set, metadata_hash=self._compute_hash_from_metadata(self.metadata))
        return [target]







class FeatureTableInput(Input):
    # PRELIMINARY ROUGH IMPLEMENTATION

    """
    Input subclass for handling feature table.
    """

    def __init__(self, dataset, preproc, samples_to_add, exclude_samples, additional_input_options ):
        """
        Initializes a FeatureTableSpecificationInput instance.

        Args:
            dataset (Dataset): The dataset associated with the feature table.
            sample_set (str): The name of the sample set used to create the feature table.
            feature_table_name (str): User-defined or default name for the feature table.
            source_wc_string (str): Wildcard string representing the original produced feature table.
        """
        super().__init__(dataset, None)
        self.additional_input_options = additional_input_options
        self.feature_table_name = additional_input_options['feature_table_name']
        self.source_wc_string = additional_input_options['source_wc_string']
        self.sample_set = additional_input_options['sample_set_name']



    def get_input_list_for_formatting(self):

        formatting_list = [
            dict(
                fs_prefix=self.dataset.fs_prefix,
                df=self.dataset.dataset_name,
                ft_name = self.feature_table_name,  
                sample_set=self.sample_set
            )
        ]
        return formatting_list

