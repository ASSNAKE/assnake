import os
import click
import glob

from assnake.core.preset_manager import PresetManager

from assnake.utils.general import read_yaml
from assnake.core.sample_set import generic_command_individual_samples, generate_result_list, generic_command_dict_of_sample_sets, prepare_sample_set_tsv_and_get_results
from assnake.core.command_builder import sample_set_construction_options, add_options

from pkg_resources import iter_entry_points

from assnake.core.snake_module import SnakeModule


class Result:
    '''
    This class is used to create InvocationCommand for CLI.

    '''
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

    def generate_cli_command(self):
        preset_options = []
        if self.preset_manager is not None:
            preset_options = self.preset_manager.gen_click_option()

        if self.input_type == 'illumina_strand_file':
            @click.command(self.name, short_help=self.description)
            @add_options(sample_set_construction_options)
            @click.option('--strand', help='Strand to profile. Default - R1', default='R1', type=click.STRING)
            @add_options(preset_options)
            @add_options(self.additional_inputs)
            @click.pass_obj
            def result_invocation(config, strand, **kwargs):
                sample_set, sample_set_name = generic_command_individual_samples(config,  **kwargs)
                config['requests'] += generate_result_list(sample_set, self.result_wc, strand=strand)

            return result_invocation
        elif self.input_type == 'illumina_sample':
        
            @click.command(self.name, short_help=self.description)
            @add_options(sample_set_construction_options)
            @add_options(preset_options)
            @add_options(self.additional_inputs)
            @click.pass_obj
            def result_invocation(config, **kwargs):
                # kwargs.update({'params': params})
                sample_set, sample_set_name = generic_command_individual_samples(
                    config,  **kwargs)
                config['requests'] += generate_result_list(
                    sample_set, self.result_wc, **kwargs)

            return result_invocation

        elif self.input_type == 'illumina_strand_file_set':
        
            @click.command(self.name, short_help=self.description)
            @add_options(sample_set_construction_options)
            @add_options(preset_options)
            @add_options(self.additional_inputs)
            @click.option('--strand', help='Strand to profile. Default - R1', default='R1', type=click.STRING)
            @click.pass_obj
            def result_invocation(config, strand, **kwargs):
                sample_sets = generic_command_dict_of_sample_sets(config,  **kwargs)
                sample_set_dir_wc = self.wc_config[self.name+'_strand_file_set_dir_wc']
                result_wc = self.wc_config[self.name + '_wc']
                res_list = prepare_sample_set_tsv_and_get_results(
                    sample_set_dir_wc, result_wc, df=kwargs['df'], sample_sets=sample_sets, strand=strand, overwrite=False)
                config['requests'] += res_list

            return result_invocation


        elif self.input_type == 'illumina_sample_set':
        
            @click.command(self.name, short_help=self.description)
            @add_options(sample_set_construction_options)
            @add_options(preset_options)
            @add_options(self.additional_inputs)
            @click.pass_obj
            def result_invocation(config, **kwargs):
                sample_sets = generic_command_dict_of_sample_sets(config,   **kwargs)

                sample_set_dir_wc = self.wc_config[self.name+'_sample_set_tsv_wc']
                result_wc = self.wc_config[self.name + '_wc']
                res_list = prepare_sample_set_tsv_and_get_results(sample_set_dir_wc, result_wc, sample_sets = sample_sets,  **kwargs)

                config['requests'] += res_list

            return result_invocation

    @classmethod
    def from_location(cls, name, result_type, location, input_type,  additional_inputs= [], description = '', with_presets = False, static_files_dir_name = 'static', invocation_command = None, preset_preparation = None):
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
                dir_in_database = os.path.join('presets', name),
                included_presets_dir = os.path.join(location, 'presets'),
                static_files_dir_name = static_files_dir_name)
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
        # Discover plugins
        discovered_plugins = SnakeModule.get_all_modules_as_dict()
        for module_name, module_class in discovered_plugins.items():
            print(module_name)
            for res in module_class.results:

                # print('\t' + str(res.invocation_command))

                if res.preset_preparation is not None:
                    print('\t' + res.name)
                    print('\t' + str(res.preset_preparation))
                    res.preset_preparation()
            # config.update({module_name:module_class.install_dir})
            # for wc_conf in module_class.wc_configs:
            #     if wc_conf is not None:
            #         wc_config.update(wc_conf)

        return discovered_plugins
