import os

class Result:

    name = '' # This name will be used as command name and in web-service
    result_wc = '' # Just ONE wc that will be populated with data and passed to snakemake, so it can understand what it needs to do
    input_type = '' # sample_file, sample, sample_set ...
    additional_inputs = '' # db, fasta_ref ...
    git_results = [] # List of files that will be added to git repo
    workflows = [] # List of files with snakemake workflow files (.smk)


    def __init__(self, name, result_wc, input_type, additional_inputs, git_results, workflows):
        self.name = name
        self.result_wc = result_wc
        self.input_type = input_type
        self.additional_inputs = additional_inputs
        self.git_results = git_results
        self.workflows = workflows

    def generate_result_list(self, input, additional_inputs = None):
        res_list = []

        return res_list