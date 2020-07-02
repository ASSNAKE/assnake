import os, click, glob

from assnake.utils.general import read_yaml
from assnake.core.sample_set import generic_command_individual_samples, generate_result_list
from assnake.core.command_builder import sample_set_construction_options, add_options

from pkg_resources import iter_entry_points 

from assnake.core.snake_module import get_all_modules_as_dict

class Result:
    '''
    This class is used to create InvocationCommand for CLI.
     
    '''
    name = '' # This name will be used as command name and in web-service
    result_wc = '' # Just ONE wc that will be populated with data and passed to snakemake, so it can understand what actions it needs to take
    input_type = '' # sample_file, sample, sample_set ...
    additional_inputs = ''  # db, fasta_ref ...
    git_results = [] # List of files that will be added to git repo
    workflows = [] # List of files with snakemake workflow files (.smk)

    invocation_command = None
    params_preparation = None # Functions that defines how to prepare and store params for this result. On installation, it checks if assnake is configured, and deploys params and necessary static files to the database. If assnake is not configured, this commands will be run during configuration.

    def __init__(self, name, workflows, result_wc, input_type, additional_inputs, invocation_command = None, wc_config = None, params_preparation = None): 
        self.name = name
        self.result_wc = result_wc
        self.input_type = input_type
        self.additional_inputs = additional_inputs
        self.workflows = workflows
        self.wc_config = wc_config
        self.params_preparation = params_preparation

        if invocation_command == None:
            if self.input_type == 'illumina_strand_file':
                @click.command(self.name, short_help=self.name)
                @add_options(sample_set_construction_options)
                @click.pass_obj
                def result_invocation(config, **kwargs):
                    kwargs.update({'strand': 'R1'})
                    sample_set, sample_set_name = generic_command_individual_samples(config,  **kwargs)
                    config['requests'] += generate_result_list(sample_set, self.result_wc, **kwargs)

                self.invocation_command = result_invocation
            elif self.input_type == 'illumina_sample':
                @click.command(self.name, short_help=self.name)
                @add_options(sample_set_construction_options)
                @click.pass_obj
                def result_invocation(config, params, **kwargs):
                    kwargs.update({'params': params})
                    sample_set, sample_set_name = generic_command_individual_samples(config,  **kwargs)
                    config['requests'] += generate_result_list(sample_set, self.result_wc, **kwargs)
        elif callable(invocation_command):
            self.invocation_command = invocation_command

    @classmethod
    def from_location(cls, name, location, input_type, additional_inputs = None, invocation_command = None, params_preparation = None):

        # Try to find all needed files for Result
        workflows = glob.glob(os.path.join(location, 'workflow*.smk'))
        wc_config = os.path.join(location, 'wc_config.yaml')

        if len(workflows) > 0:
            pass
        else:
            print('===', name)
            print('No workflow files found')

        if os.path.isfile(wc_config):
            wc_config = read_yaml(wc_config)
            result_wc = wc_config[name + '_wc']
        else:
            print('No wc config file found for result', name)
            wc_config = {}
            result_wc = ''

        x = cls(
            name = name, 
            workflows = workflows,
            result_wc = result_wc,
            input_type =  input_type,
            additional_inputs = additional_inputs,
            invocation_command = invocation_command,
            params_preparation = params_preparation,
            wc_config = wc_config
        )

        return x

    def generate_result_list(self, input, additional_inputs = None):
        res_list = []

        return res_list


def get_all_results_as_list():
    # Discover plugins
    discovered_plugins = get_all_modules_as_dict()
    for module_name, module_class in discovered_plugins.items():
        print(module_name)
        for res in module_class.results:
            
            # print('\t' + str(res.invocation_command))

            if res.params_preparation is not None:
                print('\t' + res.name)
                print('\t' + str(res.params_preparation))
                res.params_preparation()
        # config.update({module_name:module_class.install_dir})
        # for wc_conf in module_class.wc_configs:
        #     if wc_conf is not None:
        #         wc_config.update(wc_conf)
        
    return discovered_plugins