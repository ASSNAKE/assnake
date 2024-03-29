from datetime import datetime
import os
import click
from assnake.core.config import read_internal_config
from assnake.core.SampleContainerSet import SampleContainerSet
from assnake.core.Dataset import Dataset
from assnake.core.Result import Result, Step


class Pipeline:
    """
    Represents a pipeline in the Assnake system, capable of handling both preprocessing and 
    analytical steps for a given dataset. Allows dynamic configuration of processing steps 
    with specified results and presets.

    NOT FINISHED CLASS. NEEDS TO BE UPDATED WITH STEP AND CHAIN INTEGRATION
    """

    def __init__(self, name, description, preprocessing_chain=None, analytical_chain=None):
        self.name = name
        self.description = description
        self.preprocessing_chain = preprocessing_chain if preprocessing_chain is not None else {}
        self.analytical_chain = analytical_chain if analytical_chain is not None else {}
        self.steps = []



    def add_result(self, result_name, **kwargs):
        result = Result.get_result_by_name(result_name)
        if result:
            step = result.create_step(kwargs, kwargs)
            self.steps.append(step)
        return self

    def get_snakemake_targets(self):
        targets = []
        for step in self.steps:
            targets.extend(step.prepare_targets())
        return targets


    def set_dataset(self, dataset_name):
        self.dataset = Dataset(dataset_name)

    def construct_target_path(self, previous_step, result_name, preset):
        # Construct the target path based on the given parameters
        return f"{previous_step}__{result_name}_{preset}"

    def prepare_final_preprocessing_targets(self):
        targets = []
        previous_step = "raw"

        for step_number, step_details in sorted(self.preprocessing_chain.items()):
            result_name = step_details["result"]
            preset = step_details["default_preset"]

            res = Result.get_result_by_name(result_name)

            target = self.construct_target_path(previous_step, res.preprocessing_name_in_wc, preset)
            targets.append(target)
            previous_step = target

        return targets
    
    def prepare_preprocessing_chain_targets(self):
        """
        Constructs target paths for the final output files of the preprocessing chain for each sample in the dataset.

        Returns:
            List[str]: A list of target file paths for the final preprocessing outputs.
        """
        # Determine the final preprocessing step
        final_preprocessing_step = self.prepare_final_preprocessing_targets()[-1]

        sample_container_set = SampleContainerSet.create_from_kwargs(self.dataset.dataset_name, 'raw', [], [])


        # Construct target paths for each sample
        targets = []
        for container in sample_container_set.sample_containers:
            target_path = f"{container.fs_prefix}/{container.dataset_name}/reads/{final_preprocessing_step}/{container.sample_id}_R1.fastq.gz"
            targets.append(target_path)

        return targets
    


    def execute(self, config, snakemake_targets, run):
        import snakemake # Moved import here because it is slow as fucking fuck
    
        internal_config = read_internal_config()

        click.secho('-----===RUN SNAKEMAKE===-----', bg='green', fg='black')

        curr_dir = os.path.abspath(os.path.dirname(__file__))
        # click.echo('Current dir: ' + curr_dir)
        # click.echo('Number of jobs torun in parallel: ' + str(jobs))

        config['requests'] += snakemake_targets
        status = snakemake.snakemake(os.path.join(curr_dir, '../snake/snake_base.py'), 
            config = config['wc_config'],    
            targets=config['requests'], 
            printshellcmds=True,
            dryrun=not run, 
            # config = load_config_file(),
            configfiles=[internal_config['instance_config_loc']],
            drmaa_log_dir = config['config']['drmaa_log_dir'],
            

            conda_frontend='conda',
            use_conda = True,
            rerun_triggers = ['mtime'],
            conda_prefix = config['config']['conda_dir'],
            keepgoing = True,
            latency_wait = 10,
            cores=100, nodes=100)

        # print(config['requested_results']) 
        
        return status

