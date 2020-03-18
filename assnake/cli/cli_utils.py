import assnake.api.loaders
import assnake
from tabulate import tabulate
import click, os
import pandas as pd

# https://stackoverflow.com/a/40195800
options_wo_params = [
    click.option('--df','-d', help='Name of the dataset', required=True ),
    click.option('--preproc','-p', help='Preprocessing to use' ),

    click.option('--meta-column', '-c', help='Select samples based on metadata column' ),
    click.option('--column-value','-v', help='Value of metadata column by which select samples' ),

    click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING ),
    click.option('--exclude-samples','-x', 
                help='Exclude this samples from run', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING ),
]

options_w_params = [
    click.option('--df','-d', help='Name of the dataset', required=True ),
    click.option('--preproc','-p', help='Preprocessing to use' ),

    click.option('--meta-column', '-c', help='Select samples based on metadata column' ),
    click.option('--column-value','-v', help='Value of metadata column by which select samples' ),

    click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING ),
    click.option('--exclude-samples','-x', 
                help='Exclude this samples from run', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING ),
    click.option('--params', 
                help='Parameters to use', 
                default='def',
                type=click.STRING )
]


def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options


def generic_command_individual_samples(config, wc_str, df, preproc, meta_column, column_value, samples_to_add, exclude_samples, params):
    exclude_samples = [] if exclude_samples == '' else [c.strip() for c in exclude_samples.split(',')]

    samples_to_add = [] if samples_to_add == '' else [c.strip() for c in samples_to_add.split(',')]
    df_loaded = assnake.api.loaders.load_df_from_db(df)
    config['requested_dfs'] += [df_loaded['df']]

    # Now for the meta column stuff
    meta_loc = os.path.join(df_loaded['fs_prefix'], df_loaded['df'], 'mg_samples.tsv')
    if os.path.isfile(meta_loc):
        meta = pd.read_csv(meta_loc, sep = '\t')
        if meta_column is not None and column_value is not None:
            subset_by_col_value = meta.loc[meta[meta_column] == column_value]
            if len(subset_by_col_value) > 0:
                samples_to_add = list(subset_by_col_value['sample_name'])


    sample_set = assnake.SampleSet(df_loaded['fs_prefix'], df_loaded['df'], preproc, samples_to_add=samples_to_add)
    res_list = []

    if len(exclude_samples) > 0 :  
        sample_set.samples_pd = sample_set.samples_pd.loc[~sample_set.samples_pd['fs_name'].isin(exclude_samples), ]

    click.echo(tabulate(sample_set.samples_pd[['fs_name', 'reads', 'preproc']].sort_values('reads'), headers='keys', tablefmt='fancy_grid'))

    for s in sample_set.samples_pd.to_dict(orient='records'):
        preprocessing = s['preproc']
        res_list.append(wc_str.format(
            fs_prefix = s['fs_prefix'].rstrip('\/'),    
            df = s['df'],
            preproc = preprocessing,
            sample = s['fs_name'],
            params = params
        ))

    if config.get('requests', None) is None:
        config['requests'] = res_list
    else:
        config['requests'] += res_list
    




# def magic_options(func):
#     @add_options(options)
#     @click.pass_obj
#     def distill_magic(config, **kwargs):
#         kwargs['result'] = 'ssss'
#         func(config, **kwargs)

#     return distill_magic

# @click.command('test')
# @magic_options
# def tt(config, **kwargs):
#     print(kwargs)