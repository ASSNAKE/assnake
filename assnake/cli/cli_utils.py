

import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import click

@click.command('remove-human-bbmap', short_help='Count number of reads and basepairs in fastq file')

@click.option('--df','-d', help='Name of the dataset', required=True )
@click.option('--preproc','-p', help='Preprocessing to use' )

@click.option('--meta-column', '-c', help='Select samples based on metadata column' )
@click.option('--column-value','-v', help='Value of metadata column by which select samples' )

@click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING )
@click.option('--exclude-samples','-x', 
                help='Exclude this samples from run', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING )

@click.pass_obj

def generic_command_individual_samples(config, df, preproc, meta_column, meta_value, samples_to_add, exclude_samples):
    samples_to_add = [] if samples_to_add == '' else [c.strip() for c in samples_to_add.split(',')]
    exclude_samples = [] if exclude_samples == '' else [c.strip() for c in exclude_samples.split(',')]

    df = assnake.api.loaders.load_df_from_db(df)
    config['requested_dfs'] += [df['df']]
    ss = assnake.api.sample_set.SampleSet(df['fs_prefix'], df['df'], preproc, samples_to_add=samples_to_add)

    click.echo(tabulate(ss.samples_pd[['fs_name', 'reads', 'preproc']].sort_values('reads'), 
        headers='keys', tablefmt='fancy_grid'))
    res_list = []

    if len(exclude_samples) > 0 :  
        ss.samples_pd = ss.samples_pd.loc[~ss.samples_pd['fs_name'].isin(exclude_samples), ]
    for s in ss.samples_pd.to_dict(orient='records'):
        preprocessing = s['preproc']
        res_list.append( '{fs_prefix}/{df}/reads/{preproc}__rmhum_bbmap/{sample}_R1.fastq.gz'.format(
            fs_prefix = s['fs_prefix'].rstrip('\/'),
            df = s['df'],
            preproc = preprocessing,
            sample = s['fs_name']
        ))

    if config.get('requests', None) is None:
        config['requests'] = res_list
    else:
        config['requests'] += res_list
