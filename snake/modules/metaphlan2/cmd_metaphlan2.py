import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import click, glob

parameters = [p.split('/')[-1].replace('.json', '') for p in glob.glob('/data11/bio/databases/ASSNAKE/params/mp2/*.json')]

@click.command('trimmomatic', short_help='Quality based trimming')

@click.option('--df','-d', help='Name of the dataset', required=True )
@click.option('--preproc','-p', help='Preprocessing to use' )
@click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING )
@click.option('--params', help='Parameters id to use. Available parameter sets: ' + str(parameters), required=False, default = 'def')


@click.pass_obj

def cli(config, df, preproc, samples_to_add, params):
    samples_to_add = [] if samples_to_add == '' else [c.strip() for c in samples_to_add.split(',')]
    df = assnake.api.loaders.load_df_from_db(df)
    ss = assnake.api.sample_set.SampleSet(df['fs_prefix'], df['df'], preproc, samples_to_add=samples_to_add)

    click.echo(tabulate(ss.samples_pd[['fs_name', 'reads', 'preproc']].sort_values('reads'), 
        headers='keys', tablefmt='fancy_grid'))
    res_list = ss.get_locs_for_result('mp2')
    if config.get('requests', None) is None:
        config['requests'] = res_list
    else:
        config['requests'] += res_list