import os
import click
import glob
from assnake.core.config import read_assnake_instance_config
from assnake.new_core.SampleContainerSet import SampleContainerSet
from assnake.new_core.PresetManager import PresetManager

from assnake.utils.general import read_yaml
from assnake.core.command_builder import sample_set_construction_options, add_options, custom_help

from pkg_resources import iter_entry_points

from assnake.core.snake_module import SnakeModule


class Result:
    """
    A class used to define and manage a result in the Assnake framework, primarily for CLI invocation and Snakemake integration.

    This class encapsulates the information and functionalities needed to process biological sequencing data results. It includes
    mechanisms for command-line interface generation, handling various input types, formatting wildcard strings for Snakemake, and 
    managing sample sets and presets.

    Attributes:
        name (str): The name of the result, used in the CLI and web service.
        result_wc (str): A wildcard string used by Snakemake to understand required actions.
        input_type (str): The type of primary input for the result (e.g., 'sample_file', 'sample', 'sample_set').
        additional_inputs (str): Additional input parameters (e.g., database references, fasta files).
        git_results (list): A list of files to be added to a git repository.
        workflows (list): A list of Snakemake workflow files (with '.smk' extension).
        invocation_command (function): A callable that generates the CLI command for the result.
        preset_preparation (function): A function defining how to prepare and store parameters for this result.
        preset_manager (PresetManager): An instance of PresetManager to handle preset configurations.

    Methods:
        generate_cli_command: Generates a Click command for the result based on its configuration.
        format_result_wc: Formats wildcard strings for Snakemake based on sample containers and presets.
        prepare_sample_set_tsv_and_get_target_file: Prepares sample sets in TSV format and generates corresponding Snakemake targets.
        from_location (classmethod): Creates an instance of Result by automatically locating necessary files and configurations.

    Static Methods:
        get_all_results_as_list: Retrieves all Result objects from the discovered plugins.
        get_result_by_name: Finds and returns a Result object by its name.

    Usage:
        The Result class is typically instantiated with specific parameters defining a result type in the Assnake system.
        It is used to generate command-line interfaces for specific analytical results and to integrate with Snakemake workflows.
    """
    
    name = ''  # This name will be used as command name and in web-service
    result_wc = ''  # Just ONE wc that will be populated with data and passed to snakemake, so it can understand what actions it needs to take
    input_type = ''  # sample_file, sample, sample_set ...
    additional_inputs = ''  # db, fasta_ref ...
    git_results = []  # List of files that will be added to git repo
    workflows = []  # List of files with snakemake workflow files (.smk)

    invocation_command = None
    preset_preparation = None  # Functions that defines how to prepare and store params for this result. On installation, it checks if assnake is configured, and deploys params and necessary static files to the database. If assnake is not configured, this commands will be run during configuration.

    def __init__(self, name, result_type,  workflows, result_wc, input_type, additional_inputs, description = '', invocation_command=None, wc_config=None, preset_preparation=None, preset_manager = None):
        self.name = name
        self.result_type = result_type
        self.result_wc = result_wc
        self.sample_set_dir_wc = wc_config.get(f'{name}_sample_set_dir_wc', None)
        self.input_type = input_type
        self.additional_inputs = additional_inputs
        self.description = description
        self.workflows = workflows
        self.wc_config = wc_config
        self.preset_preparation = preset_preparation
        self.preset_manager = preset_manager

        if invocation_command == None:
            self.invocation_command = self.generate_cli_command()
        elif callable(invocation_command):
            self.invocation_command = invocation_command()

    def handle_illumina_strand_file(self, config, kwargs):
        sample_container_set = SampleContainerSet.create_from_kwargs(**{key: kwargs[key] for key in ['df', 'preproc', 'samples_to_add', 'exclude_samples']})
        config['requests'] += self.format_result_wc(sample_container_set, strand=kwargs['strand'])

    def handle_illumina_sample(self, config, kwargs):
        sample_container_set = SampleContainerSet.create_from_kwargs(**{key: kwargs[key] for key in ['df', 'preproc', 'samples_to_add', 'exclude_samples']})
        config['requests'] += self.format_result_wc(sample_container_set, kwargs['preset'])

    def handle_illumina_strand_file_set(self, config, kwargs):
        sample_container_set = SampleContainerSet.create_from_kwargs(**{key: kwargs[key] for key in ['df', 'preproc', 'samples_to_add', 'exclude_samples']})
        formatted_target = self.prepare_sample_set_tsv_and_get_target_file(sample_container_set, kwargs['sample_set'], overwrite=kwargs.get('overwrite', False), **kwargs)
        config['requests'] += [formatted_target]

    def handle_illumina_sample_set(self, config, kwargs):
        sample_container_set = SampleContainerSet.create_from_kwargs(**{key: kwargs[key] for key in ['df', 'preproc', 'samples_to_add', 'exclude_samples']})
        formatted_target = self.prepare_sample_set_tsv_and_get_target_file(sample_container_set, kwargs['sample_set'], overwrite=kwargs.get('overwrite', False), **kwargs)
        config['requests'] += [formatted_target]

    def generate_cli_command(self):
        # Common setup for all command types
        preset_options = self.preset_manager.gen_click_option() if self.preset_manager is not None else []
        strand_option = [click.option('--strand', help='Strand to profile. Default - R1', default='R1', type=click.STRING)] if self.input_type in ['illumina_strand_file', 'illumina_strand_file_set'] else []

        # Unified command decorator
        @click.command(self.name, short_help=self.description, add_help_option=False)
        @add_options(sample_set_construction_options)
        @add_options(preset_options)
        @add_options(self.additional_inputs)
        @add_options(strand_option)
        @click.pass_obj
        def result_invocation(config, **kwargs):
            if not kwargs.get('df'):
                custom_help(ctx=click.get_current_context(), param=None, value=True)

            if 'preset' in kwargs and self.preset_manager is not None:
                preset = self.preset_manager.find_preset_by_name_and_dataset(kwargs['preset'], kwargs['df'])
                kwargs['preset'] = preset.name_wo_ext
                    
            if self.input_type == 'illumina_strand_file':
                self.handle_illumina_strand_file(config, kwargs)
            elif self.input_type == 'illumina_sample':
                self.handle_illumina_sample(config, kwargs)
            elif self.input_type == 'illumina_strand_file_set':
                self.handle_illumina_strand_file_set(config, kwargs)
            elif self.input_type == 'illumina_sample_set':
                self.handle_illumina_sample_set(config, kwargs)
            else:
                click.secho(f'Unsupported input type: {self.input_type}', fg='red')
                exit()
            
        return result_invocation

    def format_wildcard_string(self, container, preset, additional_inputs):
        # Shared logic for formatting wildcard strings
        formatting_dict = {
            'fs_prefix': container.fs_prefix,
            'df': container.dataset_name,
            'preproc': container.preprocessing,
            'df_sample': container.sample_id,
            'preset': preset
        }
        formatting_dict.update(additional_inputs)
        return self.result_wc.format(**formatting_dict)

    def format_result_wc(self, sample_container_set, preset='default', **additional_inputs):
        print(additional_inputs)
        return [self.format_wildcard_string(container, preset, additional_inputs) for container in sample_container_set.sample_containers]

    @property
    def preprocessing_name_in_wc(self):
        if self.result_type == 'preprocessing':
            preprocessing_name = self.result_wc.split('/')[3]
            preprocessing_name = preprocessing_name.replace('{preproc}__', '').replace('_{preset}', '')
            return preprocessing_name
        else:
            return 0

    
    def prepare_sample_set_tsv_and_get_target_file(self, sample_container_set, sample_set_name, overwrite=False, **kwargs):
        """
        Prepares sample sets in TSV format and generates a list of results for Snakemake.

        Args:
            sample_container_set (SampleContainerSet): The SampleContainerSet to be saved and processed.
            sample_set_name (str): Name of the sample set.
            overwrite (bool): Flag to indicate whether existing files should be overwritten.
            **kwargs: Additional keyword arguments to be included in the result wildcard strings.

        Returns:
            list: A list of formatted result wildcard strings for Snakemake.
        """
        # Construct the directory path for the sample set
        sample_set_dir = self.sample_set_dir_wc.format(
            fs_prefix=sample_container_set.dataset.fs_prefix,
            df=sample_container_set.dataset.dataset_name,
            sample_set=sample_set_name,
            **kwargs
        )
        sample_set_loc = os.path.join(sample_set_dir, 'sample_set.tsv')

        # Create the directory if it does not exist
        if not os.path.exists(sample_set_dir):
            os.makedirs(sample_set_dir, exist_ok=True)

        # Check if the sample set file exists
        if not os.path.isfile(sample_set_loc) or overwrite:
            if os.path.isfile(sample_set_loc):
                click.secho('Overwriting...', fg='yellow')
            sample_container_set.to_tsv(sample_set_loc, overwrite) 
            
        else:
            click.secho('Sample set with this name already exists!', fg='red')

        # Generate and return the list of formatted result wildcard strings
        formatting_dict = {
                'fs_prefix': sample_container_set.dataset.fs_prefix,
                'df': sample_container_set.dataset.dataset_name,
                'sample_set': sample_set_name
            }
        formatting_dict.update(kwargs)

        return self.result_wc.format(**formatting_dict)


    @classmethod
    def from_location(cls, name, result_type, location, input_type,  additional_inputs= [], description = '', with_presets = False, static_files_dir_name = 'static', invocation_command = None, preset_preparation = None, preset_file_format = 'json'):
        '''
        This method is a wrapper for creating an instance of Result class by grabbing the snakefiles and wc_config automatically 
        from user-provided `location`. How can we make grabbing the location automatic? 

        :param name:       Name of the Result, will also appear as a command name.
        :param location:   Absolute location where files are stored. We need it in order to grab .smk and wc_config and all other static files.
        :param input_type: Type of primary input for this result. Currently we use illumina_sample, illumina_sample_set, illumina_strand_file, illumina_strand_file_set 

        :param additional_inputs:  Any additional inputs for the Result. Usually it would be Database, Reference and their mixes.
        :param invocation_command: Users can pass custom invocation commands that will be used in CLI instead of auto-generated ones. Useful for custom logic in CLI.
        :param preset_preparation: Functions that tells the Result how to save it's presets to database. 
        :param preset

        :param description: Description of the Result, will also appear as a short help in CLI.
        '''
        # Find all smk files for Result. This is snakemake workflows.
        workflows = glob.glob(os.path.join(location, 'workflow*.smk'))
        # Find all wc_config files.
        wc_config = os.path.join(location, 'wc_config.yaml')
        # Find default config files
        default_config = os.path.join(location, 'default_config.yaml')

        instance_config = read_assnake_instance_config()


        if len(workflows) == 0:  # No workflow files found
            print('===', name)
            print('No workflow files found')

        # Check for config file with wildcards and load.
        if os.path.isfile(wc_config):
            wc_config = read_yaml(wc_config)
            result_wc = wc_config[name + '_wc']
        else:
            # print('No wc config file found for result', name)
            wc_config = {}
            result_wc = ''

        if os.path.isfile(default_config):
            default_config = read_yaml(default_config)
        else:
            # print('No default_config file found for result', name)
            default_config = None

        if with_presets:
            preset_manager = PresetManager(
                dir_in_database=os.path.join(instance_config['assnake_db'], 'presets', name),
                included_presets_dir=os.path.join(location, 'presets'),
                static_files_dir_name=static_files_dir_name,
                preset_file_format=preset_file_format,
                module_name='NOT READY NEED TO UPDATE RESULT CREATION IN MODULES',  # Add the module name
                result_name=name  # Add the result name
            )
        else:
            preset_manager = None

        x = cls(
            name = name,
            workflows= workflows,
            result_type = result_type,
            result_wc= result_wc,
            input_type=  input_type,
            additional_inputs= additional_inputs,
            invocation_command= invocation_command,
            preset_preparation= preset_preparation,
            wc_config= wc_config,
            preset_manager = preset_manager,
            description=description,
        )

        return x
    

    

    @staticmethod
    def get_all_results_as_list():
        """
        Retrieves a list of all Result objects from the discovered plugins.

        This method iterates through all modules discovered by SnakeModule, 
        aggregates their results, and returns a list of these Result objects.

        Returns:
            List[Result]: A list containing all Result objects from the discovered plugins.
        """
        results = []
        discovered_plugins = SnakeModule.get_all_modules_as_dict()  # Discover all plugins

        # Iterate through each discovered plugin and append their results to the list
        for module_name, module_class in discovered_plugins.items():
            print('module_name', module_name)
            for res in module_class.results:
                print('\t' + res.name)  # Logging the name of each result for debugging
                results.append(res)

        return results
    
    @staticmethod
    def get_result_by_name(result_name):
        """
        Finds and returns a Result object by its name from the discovered plugins.

        This method searches through all modules discovered by SnakeModule to find
        a Result object matching the provided name.

        Args:
            result_name (str): The name of the Result object to search for.

        Returns:
            Result or None: The Result object with the matching name, if found; otherwise, None.
        """
        discovered_plugins = SnakeModule.get_all_modules_as_dict()  # Discover all plugins

        # Iterate through each discovered plugin to find the Result with the specified name
        for module_name, module_class in discovered_plugins.items():
            for res in module_class.results:
                if res.name == result_name:
                    return res  # Return the found Result object

        return None  # Return None if no matching Result is found
    

    def __repr__(self):
        return self.name
    
    def __str__(self) -> str:
        attributes = vars(self)
        final = ''
        for attr_name, attr_value in attributes.items():
            final += f"{attr_name}: {attr_value}"
            final += '\n'
