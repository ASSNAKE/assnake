import os, click, glob
from assnake.cli.cli_utils import generic_command_individual_samples, generate_result_list, sample_set_construction_options, add_options
from assnake.utils import read_yaml

class Result:

    name = '' # This name will be used as command name and in web-service
    result_wc = '' # Just ONE wc that will be populated with data and passed to snakemake, so it can understand what it needs to do
    input_type = '' # sample_file, sample, sample_set ...
    additional_inputs = '' # db, fasta_ref ...
    git_results = [] # List of files that will be added to git repo
    workflows = [] # List of files with snakemake workflow files (.smk)

    imvocation_command = None

    def __init__(self, name, workflows, result_wc, input_type, additional_inputs, invocation_command = None, wc_config = None): 
        self.name = name
        self.result_wc = result_wc
        self.input_type = input_type
        self.additional_inputs = additional_inputs
        self.workflows = workflows
        self.wc_config = wc_config

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
                @click.command(self.name + '_fff', short_help=self.name)
                @add_options(sample_set_construction_options)
                @click.pass_obj
                def trimmomatic_invocation(config, params, **kwargs):
                    kwargs.update({'params': params})
                    sample_set, sample_set_name = generic_command_individual_samples(config,  **kwargs)
                    config['requests'] += generate_result_list(sample_set, self.result_wc, **kwargs)
        elif callable(invocation_command):
            self.invocation_command = invocation_command

    @classmethod
    def from_location(cls, name, location, input_type, additional_inputs = None, invocation_command = None):

        # Try to find all needed files for Result
        workflows = glob.glob(os.path.join(location, 'workflow*.smk'))
        wc_config = os.path.join(location, 'wc_config.yaml')

        if len(workflows) > 0:
            [print(w) for w in workflows]
        else:
            print('===', name)
            print('No workflow files found')

        if os.path.isfile(wc_config):
            wc_config = read_yaml(wc_config)
        else:
            print('No config file found')
            exit()

        x = cls(
            name = name, 
            workflows = workflows,
            result_wc = wc_config[name + '_wc'],
            input_type =  input_type,
            additional_inputs = additional_inputs,
            invocation_command = invocation_command,
            wc_config=wc_config
        )

        return x

    def generate_result_list(self, input, additional_inputs = None):
        res_list = []

        return res_list

    