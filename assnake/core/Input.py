from datetime import datetime
import os
import yaml
import hashlib

import re

from assnake.core.SampleContainerSet import SampleContainerSet
from assnake.core.Dataset import Dataset

class Input:
    def __init__(self, dataset, preprocessing=None, samples_to_add=None, exclude_samples=None, additional_input_options=None):
        self.dataset = dataset
        self.dataset = Dataset(self.dataset)

        self.preprocessing = preprocessing
        self.samples_to_add = samples_to_add
        self.exclude_samples = exclude_samples
        self.additional_input_options = additional_input_options or {}
        if preprocessing:
            self.sample_container_set = SampleContainerSet(self.dataset, preprocessing, samples_to_add, exclude_samples)
        else:
            self.sample_container_set = None


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
    
    
class FeatureTableSpecificationInput(Input):
    """
    Input subclass for handling feature table specifications.
    """

    def __init__(self, dataset, additional_input_options):
        super().__init__(dataset, None)
        self.feature_table_name = additional_input_options.get('feature_table_name')
        self.source_wc_string = additional_input_options.get('source_wc_string')
        self.sample_set = additional_input_options.get('sample_set')
        self.parsed_presets = additional_input_options.get('parsed_presets')
        # Dynamically parse presets and other info from source_wc_string
        self.parsed_metadata = self.parse_dynamic_metadata(self.source_wc_string, additional_input_options)
        
        # Save and retrieve metadata
        self.metadata = self._create_metadata()
        self.feature_table_path = self._get_feature_table_path()

    def parse_dynamic_metadata(self, wc_string, options):
        # Initialize an empty metadata dictionary
        metadata = {}

        # Ensure wc_string is a valid string
        if wc_string and isinstance(wc_string, str):
            # Extract all wildcards from the source_wc_string
            wildcards = re.findall(r'\{(\w+)\}', wc_string)
            metadata = {wc: options.get(wc) for wc in wildcards}

        # Check if parsed_presets exists and update metadata if necessary
        if hasattr(self, 'parsed_presets') and self.parsed_presets:
            metadata.update(self.parsed_presets)

        return metadata

    def _create_metadata(self):
        metadata_filename = "metadata.yaml"

        metadata_filepath = os.path.join(self.dataset.full_path, 'feature_tables', self.sample_set, self.feature_table_name, metadata_filename)
        os.makedirs(os.path.dirname(metadata_filepath), exist_ok=True)
        with open(metadata_filepath, 'w') as metadata_file:
            yaml.dump(self.parsed_metadata, metadata_file)
        return self.parsed_metadata

    def _get_feature_table_path(self):
        metadata_hash = self._compute_hash_from_metadata(self.metadata)
        feature_table_dir = os.path.join(self.dataset.full_path, 'feature_tables', self.sample_set, self.feature_table_name)
        feature_table_path = os.path.join(feature_table_dir, f"{metadata_hash}.rds")
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



class FeatureTableInput(Input):
    # PRELIMINARY ROUGH IMPLEMENTATION

    """
    Input subclass for handling feature table.
    """

    def __init__(self, dataset, additional_input_options ):
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
        self.sample_set = additional_input_options['sample_set']
        self.parsed_presets = additional_input_options.get('parsed_presets')



    def get_input_list_for_formatting(self):

        formatting_list = [
            dict(
                fs_prefix=self.dataset.fs_prefix,
                df=self.dataset.dataset_name,
                ft_name = self.feature_table_name,  
                sample_set=self.sample_set,
                **self.additional_input_options['parsed_presets']
            )
        ]
        return formatting_list


class PhyloseqInput(Input):
    def __init__(self, dataset, additional_input_options ):
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
        self.sample_set = additional_input_options['sample_set']
        self.parsed_presets = additional_input_options.get('parsed_presets')
        self.filter_chain = additional_input_options.get('filter_chain')
        print(self.filter_chain)

        print(self.source_wc_string)
        if self.filter_chain == '' or self.filter_chain is None:
            self.last_step = 0
        else:
            self.last_step = int(self.filter_chain.split('/')[-1].split('_')[0].replace('step', ''))
        


    def get_input_list_for_formatting(self):

        step_num_dict = {}
        if '{step_num}' in self.source_wc_string:
            step_num_dict = {'step_num': self.last_step + 1}
        formatting_list = [
            dict(
                fs_prefix=self.dataset.fs_prefix,
                df=self.dataset.dataset_name,
                ft_name = self.feature_table_name,  
                sample_set=self.sample_set,
                filter_chain = self.filter_chain,
                **step_num_dict,
                **self.additional_input_options['parsed_presets']
            )
        ]
        return formatting_list
    



class InputFactory:
    def create_input(self, input_type, **kwargs):
        # Define expected arguments for each input type
        expected_args_for_illumina = {'dataset', 'preprocessing', 'samples_to_add', 'exclude_samples', 'additional_input_options'}
        expected_args_for_illumina_set = {'dataset', 'preprocessing', 'samples_to_add', 'exclude_samples', 'additional_input_options', 'sample_set_name'}

        expected_args_for_ft_specification = {'dataset', 'sample_set', 'feature_table_name', 'source_wc_string', 'additional_input_options', 'parsed_presets'}
        expected_args_for_ft_specification2 = {'dataset', 'additional_input_options'}

        # Ensure additional_input_options is always provided
        kwargs.setdefault('additional_input_options', {})

        # Define the expected arguments for FeatureTableSpecificationInput


        # Filter kwargs for each input type
        if input_type == 'illumina_strand_file' or input_type == 'illumina_sample':
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in expected_args_for_illumina}
            return IlluminaSampleInput(**filtered_kwargs)

        elif input_type == 'illumina_strand_file_set' or input_type == 'illumina_sample_set':
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in expected_args_for_illumina_set}
            return IlluminaSampleSetInput(**filtered_kwargs)

        if input_type == 'feature_table_specification':
            # Filter out only the expected arguments
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in expected_args_for_ft_specification2}
            # Ensure additional_input_options is a dictionary
            filtered_kwargs.setdefault('additional_input_options', {})

            additional_kwargs = {k: v for k, v in kwargs.items() if k in expected_args_for_ft_specification}

            filtered_kwargs['additional_input_options'] = additional_kwargs

            # Create the input instance
            return FeatureTableSpecificationInput(**filtered_kwargs)

        elif input_type == 'feature_table':

            filtered_kwargs = {k: v for k, v in kwargs.items() if k in expected_args_for_ft_specification2}
            additional_kwargs = {k: v for k, v in kwargs.items() if k in expected_args_for_ft_specification}
            filtered_kwargs['additional_input_options'] = additional_kwargs

            return FeatureTableInput(**filtered_kwargs)
        elif input_type == 'phyloseq':

            filtered_kwargs = {k: v for k, v in kwargs.items() if k in {'dataset', 'additional_input_options'}}
            additional_kwargs = {k: v for k, v in kwargs.items() if k in {'filter_chain', 'sample_set', 'feature_table_name', 'source_wc_string', 'additional_input_options', 'parsed_presets'}}
            filtered_kwargs['additional_input_options'] = additional_kwargs

            return PhyloseqInput(**filtered_kwargs)

        else:
            raise ValueError(f"Unsupported input type: {input_type}")


    def get_input_config(self, input_type, wc_config, name):
        # Dictionary to hold configurations for each input type
        config = {}

        if input_type == 'feature_table_specification':
            config = {
                'ft_meta_dir_wc': wc_config.get(f'{name}_ft_meta_dir_wc', None),
                'source_wc_string': wc_config.get(f'{name}_source_wc', None),
                'additional_inputs': {
                    'feature-table-name': 'Name to use for the requested feature table',
                    'sample-set': 'Name of the sample set to use. (Hash included)',
                },
                'parsed_presets': self._parse_wc_string_for_presets(wc_config.get(f'{name}_source_wc', None))
            }
        if input_type == 'feature_table':
            config = {
                'ft_meta_dir_wc': wc_config.get(f'{name}_ft_meta_dir_wc', None),
                'source_wc_string': wc_config.get(f'{name}_source_wc', None),
                'additional_inputs': {
                    'feature-table-name': 'Name to use for the requested feature table',
                    'sample-set': 'Name of the sample set to use. (Hash included)',
                },
                'parsed_presets': self._parse_wc_string_for_presets(wc_config.get(f'{name}_source_wc', None))
            }
        if input_type == 'phyloseq':
            config = {
                'ft_meta_dir_wc': wc_config.get(f'{name}_ft_meta_dir_wc', None),
                'source_wc_string': wc_config.get(f'{name}_source_wc', None),
                'additional_inputs': {
                    'feature-table-name': 'Name to use for the requested feature table',
                    'sample-set': 'Name of the sample set to use. (Hash included)',
                    'filter_chain': 'Postprocessing chain'
                },
                'parsed_presets': self._parse_wc_string_for_presets(wc_config.get(f'{name}_source_wc', None))
            }
        elif input_type == 'illumina_sample':
            # Example configuration for IlluminaSampleInput
            config = {
                'additional_inputs': {
                    # Other specific configurations for IlluminaSampleInput
                }
            }
        elif input_type == 'illumina_sample_set':
            # Example configuration for IlluminaSampleSetInput
            config = {
                'additional_inputs': {
                    'sample-set-name': 'Name of the sample set to process',
                    # Other specific configurations for IlluminaSampleSetInput
                }
            }
        # Additional input types and their configurations...

        return config
    
    def _parse_wc_string_for_presets(self, wc_string):
        # PRELIMINARY ROUGH IMPLEMENTATION

        """
        Parses the wildcard string to extract preset wildcards.

        Args:
            wc_string (str): The wildcard string to parse.

        Returns:
            dict: A dictionary where keys are the names of the presets and the values are None (to be filled later).
        """

        if wc_string is None:
            return {}
        # Regular expression pattern to match preset wildcards
        # This pattern matches anything that looks like {something_preset}
        preset_pattern = r'\{([a-zA-Z0-9_]+_preset)\}'

        # Find all matches in the wildcard string
        preset_matches = re.findall(preset_pattern, wc_string)

        try:
            preset_matches = re.findall(preset_pattern, wc_string)
            # Create a dictionary from the matches with None as default values
            presets = {preset: None for preset in preset_matches}

            return presets
        except TypeError as e:
            # Handle the exception here (e.g., print an error message, set a default value, or log the error)
            print(f"Error while parsing WC string: {e}")
            return {}
