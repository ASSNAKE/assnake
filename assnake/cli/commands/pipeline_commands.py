import click
from assnake.core.Pipeline import Pipeline

@click.command("pipeline")
@click.pass_obj
@click.option('--pipeline-script', '-s', type=click.Path(exists=True))
@click.option('--run', is_flag=True, help='Run the pipeline')
def run_pipeline(config, pipeline_script, run):
    import runpy

    # Execute the Python script and retrieve the pipeline object
    pipeline_context = runpy.run_path(pipeline_script)
    pipeline = pipeline_context['pipeline']

    # Retrieve Snakemake targets from the pipeline
    snakemake_targets = pipeline.get_snakemake_targets()

    print(snakemake_targets)

    # Run Snakemake with the targets
    pipeline.execute(config, snakemake_targets, run)



@click.command(name='test')
@click.pass_obj

@click.option('--dataset','-d', help='Name of the dataset')
@click.option('--run', is_flag=True, help='Run the pipeline')
def covid_rnf(config, dataset, run):
    my_pipeline = Pipeline(
    name="COVID RNF Pipeline", 
    description="From reads to DADA2 seqtab and tax_table in one command!",
        preprocessing_chain={
            1: {'result': 'cutadapt', 'default_preset': "RMv3v4primers"},
            2: {'result': 'dada2-filter-and-trim', 'default_preset': "strict"}
        },
        analytical_chain = {
            1: {'result':              'dada2-derep-infer-merge', 
                'sample_set_name': '18-Nov-2023_6e0ab0a3',
                'learn_errors_preset': 'def.1ac94255',
                'core_algo_preset':    'def',
                'merged_preset':       'def'},
            2: {'result':              'dada2-remove-chimeras', 
                'sample_set_name': '18-Nov-2023_6e0ab0a3',
                
                'learn_errors_preset': 'def.1ac94255',
                'core_algo_preset':    'def',
                'merged_preset':       'def',
                'nonchim_preset':      'def'},



            3: {'result':              'dada2-export-asv-table', 
                'sample_set_name': '18-Nov-2023_6e0ab0a3',
                'learn_errors_preset': 'def.1ac94255',
                'core_algo_preset':    'def',
                'merged_preset':       'def',
                'nonchim_preset':      'def',
                'ft_name': 'testft',
                'hash':    '???'},



            # 4: {'result':              'dada2-assign-taxa', 
            #     'dada2preset':         'learn_errors_def.1ac94255__core_algo__def__merged__def__nonchim__def'},
        }
    )

    my_pipeline.set_dataset(dataset)

    targets = my_pipeline.prepare_analytical_chain_targets()

    my_pipeline.execute(config, run)