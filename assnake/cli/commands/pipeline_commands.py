import click
from assnake.core.Pipeline import Pipeline

@click.command(name='test')
@click.pass_obj
@click.option('--dataset','-d', help='Name of the dataset')
@click.option('--run', is_flag=True, help='Run the pipeline')
def pipeline_testing(config, dataset, run):
    my_pipeline = Pipeline(
    name="DADA2 PE Pipeline", 
    description="From reads to DADA2 seqtab and tax_table in one command!",
        preprocessing_chain={
            1: {'result': 'cutadapt', 'default_preset': "RMv3v4primers"},
            2: {'result': 'dada2-filter-and-trim', 'default_preset': "strict"}
        },
        analytical_chain = {
            1: {'result': 'dada2-derep-infer-merge', 
                'learn_errors_preset': 'def.1ac94255',
                'core_algo_preset': 'def',
                'merged_preset': 'def'},
            2: {'result': 'dada2-remove-chimeras', 
                'learn_errors_preset': 'def.1ac94255',
                'core_algo_preset': 'def',
                'merged_preset': 'def',
                'nonchim_preset': 'def'},
            3: {'result': 'dada2-export-asv-table', 
                'learn_errors_preset': 'def.1ac94255',
                'core_algo_preset': 'def',
                'merged_preset': 'def',
                'nonchim_preset': 'def'},
            4: {'result': 'dada2-assign-taxa', 
                'dada2preset': 'learn_errors_def.1ac94255__core_algo__def__merged__def__nonchim__def'},
        }
    )

    my_pipeline.set_dataset(dataset)

    targets = my_pipeline.prepare_analytical_chain_targets()

    my_pipeline.execute(config, run)