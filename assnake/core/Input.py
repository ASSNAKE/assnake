from datetime import datetime
import os 
from assnake.core.SampleContainerSet import SampleContainerSet

class Input:
    """
    Base class for different types of inputs.
    """
    def __init__(self, dataset, preprocessing, samples_to_add, exclude_samples, additional_input_options):
        self.dataset = dataset
        self.preprocessing = preprocessing
        self.samples_to_add = samples_to_add
        self.exclude_samples = exclude_samples
        self.additional_input_options = additional_input_options
        self.sample_container_set = SampleContainerSet(dataset, preprocessing, samples_to_add, exclude_samples)

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