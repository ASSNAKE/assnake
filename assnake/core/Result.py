import os
import sys
import click
import glob
from assnake.core.config import read_assnake_instance_config
from assnake.core.SampleContainerSet import SampleContainerSet
from assnake.core.Dataset import Dataset
from assnake.core.PresetManager import PresetManager

from assnake.utils.general import read_yaml
from assnake.cli.command_builder import sample_set_construction_options, add_options, custom_help

from pkg_resources import iter_entry_points

from assnake.core.snake_module import SnakeModule
from assnake.core.Input import FeatureTableInput, IlluminaSampleInput, IlluminaSampleSetInput, FeatureTableSpecificationInput, InputFactory

from datetime import datetime
import re


class Step:
    def __init__(self, result, input_instance, other_specifications):
        """
        Initializes a Step instance.

        Args:
            result (Result): The result object associated with this step.
            input_instance (Input): The input object providing data for this step.
            other_specifications
        """
        self.result = result
        self.input_instance = input_instance
        self.other_specifications = other_specifications

    def prepare_targets(self):
        """
        Prepares the target paths for Snakemake based on the result and input.

        Returns:
            List[str]: A list of target paths for Snakemake.
        """
        return self.format_target_paths(self.input_instance, self.other_specifications)
    
    def format_target_paths(self, input_instance, formatting_dict):
        target_paths = []
        input_dicts = input_instance.get_input_list_for_formatting()

        for input_dict in input_dicts:
            input_dict.update(formatting_dict)
            target_paths.append(self.result.result_wc.format(**{**input_dict}).replace('//', '/'))
            # target_paths.append(self.result.result_wc.format(**{**input_dict}))
        return target_paths


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

    def __init__(self, name, result_type,  workflows, result_wc, input_type, additional_inputs={}, description = '', invocation_command=None, wc_config=None, preset_preparation=None, preset_manager = None):
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

        self.input_factory = InputFactory()
        self.input_config = self.input_factory.get_input_config(input_type, wc_config, name)


        self.parsed_presets = {}

        self.sample_set_construction_options = sample_set_construction_options

        # Update the result object with the configurations from the input factory
        for key, value in self.input_config.items():
            setattr(self, key, value)

        if invocation_command == None:
            self.invocation_command = self.generate_cli_command()
        elif callable(invocation_command):
            self.invocation_command = invocation_command()


    def _add_dynamic_cli_options(self, parsed_preset):
        # Dynamically add CLI options for each parsed preset wildcard
        for preset_name in self.parsed_presets:
            self.additional_inputs[preset_name] = f'{preset_name} value'

    def create_step(self, config, kwargs):
        """
        Creates a Step instance based on provided command line arguments.
        """
        # Handle preset selection and parameter override

        formatting_dict = self._handle_preset_and_params(kwargs)

        # Create input instance
        input_instance = self.input_factory.create_input(self.input_type, **{**kwargs, **self.input_config, 'sample_set_path': self.sample_set_dir_wc})

        # Create and return a Step instance
        return Step(result=self, input_instance=input_instance, other_specifications=formatting_dict)

    def _handle_preset_and_params(self, kwargs):
        """
        Handles preset file loading and parameter overrides from CLI arguments.
        """
        formatting_dict = {}

        # Habdle case when parameters are embeded into wildcard path.
        # Such wildcards have _preset suffix
        param_values_from_wc = {}
        if self.parsed_presets:
            param_values_from_wc = {param: kwargs.get(param) for param in self.parsed_presets}
        formatting_dict.update(param_values_from_wc)

        # Try to get preset from kwargs. 
        # It can only be there if result have with_presets = True
        preset_name = kwargs.get('preset', None)
        if self.preset_manager:
            # Name provided from cli, actually it should always be provided because we have default (timestamp)
            if preset_name:
                preset = self.preset_manager.find_preset_by_name_and_dataset(preset_name.split('.')[0], kwargs.get('dataset', None))
                if preset is not None:
                    formatting_dict['preset'] = preset.name_wo_ext
                else:
                    # IF WE HAVE SCHEMA
                    if self.preset_manager.schema:
                        preset_params = self.preset_manager.get_default_preset_contents()
                        for key, value in kwargs.items():
                            if key in self.preset_manager.schema.keys() and value is not None:
                                preset_params[key] = value

                        # Save new preset if parameters were overridden and a name is given
                        if preset_name and self.preset_manager:
                            self.preset_manager.save_new_preset(preset_name, preset_params)
                            preset = self.preset_manager.find_preset_by_name_and_dataset(preset_name, kwargs.get('dataset', None))
                            formatting_dict['preset'] = preset.name_wo_ext

                        # Include preset parameters in the formatting dictionary
                        # formatting_dict.update(preset_params)

        # Add any other special handling or additional parameters
        if 'strand' in kwargs:
            formatting_dict['strand'] = kwargs['strand']
        if 'sample_set_name' in kwargs:
            formatting_dict['sample_set_name'] = kwargs['sample_set_name']

        return formatting_dict
    
    def generate_cli_command(self):
        """
        Generates a Click command for the result based on its configuration.

        This method constructs a command-line interface (CLI) command using the Click library,
        allowing users to interact with the result through a command-line terminal.
        """


        if 'parsed_presets' in self.input_config:
            self._add_dynamic_cli_options(self.input_config['parsed_presets'])
        
        # Common setup for all command types
        preset_options = self.preset_manager.gen_click_option() if self.preset_manager is not None else []
        schema_options = self.preset_manager.generate_click_options_from_schema() if self.preset_manager is not None and self.preset_manager.schema is not None else []



        additional_input_options = [click.option(f'--{input_name}', help=input_description, type=click.STRING) for input_name, input_description in self.additional_inputs.items()]
        strand_option = [click.option('--strand', help='Strand to profile. Default - R1', default='R1', type=click.STRING)] if self.input_type in ['illumina_strand_file', 'illumina_strand_file_set'] else []


        # Unified command decorator
        @click.command(name=self.name, short_help=self.description, help=self.description)
        @add_options(self.sample_set_construction_options)  # Add common options for constructing a sample set

        @add_options(preset_options)                   # Add options for selecting presets
        @add_options(schema_options)  # Add options generated from the schema

        @add_options(additional_input_options)         # Add any additional input options defined for the result
        @add_options(strand_option)
        @click.pass_obj
        def result_invocation(config, **kwargs):
            """
            The actual command that will be invoked from the CLI.

            Args:
                config: The Assnake configuration object, typically passed automatically by Click.
                **kwargs: Keyword arguments capturing all the command-line options provided by the user.
            """

            # Validate dataset option
            if not kwargs.get('dataset'):
                custom_help(ctx=click.get_current_context(), param=None, value=True)


            # Create a Step instance
            step = self.create_step(config, kwargs)
            # Prepare targets using the Step instance
            snakemake_targets = step.prepare_targets()

            print('snakemake_targets', snakemake_targets)

            # Add generated targets to the configuration for Snakemake execution
            config['requests'].extend(snakemake_targets)

        return result_invocation
    

    @property
    def preprocessing_name_in_wc(self):
        if self.result_type == 'preprocessing':
            preprocessing_name = self.result_wc.split('/')[3]
            preprocessing_name = preprocessing_name.replace('{preproc}__', '').replace('_{preset}', '')
            return preprocessing_name
        else:
            return 0


    @classmethod
    def from_location(cls, name, result_type, location, input_type,  additional_inputs={}, description = '', with_presets = False, static_files_dir_name = 'static', invocation_command = None, preset_preparation = None, preset_file_format = 'json'):
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


        instance_config = read_assnake_instance_config()
        if instance_config is None:
            raise RuntimeError("Assnake instance configuration not found. Please run 'assnake config init' to configure the instance.")

        # Find all smk files for Result. This is snakemake workflows.
        workflows = glob.glob(os.path.join(location, 'workflow*.smk'))
        # Find all wc_config files.
        wc_config = os.path.join(location, 'wc_config.yaml')
        # Find default config files
        default_config = os.path.join(location, 'default_config.yaml')


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


        schema_file =  os.path.join(location, 'param_schema.yaml')

        if with_presets:
            preset_manager = PresetManager(
                dir_in_database=os.path.join(instance_config['assnake_db'], 'presets', name),
                included_presets_dir=os.path.join(location, 'presets'),
                static_files_dir_name=static_files_dir_name,
                preset_file_format=preset_file_format,
                module_name='NOT READY NEED TO UPDATE RESULT CREATION IN MODULES',
                result_name=name,
                schema_file=schema_file if os.path.isfile(schema_file) else None  # Pass the schema file path if it exists
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
