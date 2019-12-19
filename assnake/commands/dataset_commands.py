import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import snakemake

@click.command(name='list')
def df_list():
    """List datasets in database"""
    dfs = assnake.api.loaders.load_dfs_from_db('')
    if len(list(dfs.keys())) == 0:
        click.echo('No datasets in your system yet!\nYou can create one by running\n' + 
            click.style('  assnake dataset create  ', bg='blue', fg='white', bold=True))

    for df in dfs.values():
        df_name = df['df']
        click.echo(click.style(''*2 + df_name + ' '*2, fg='green', bold=True))
        click.echo('  Filesystem prefix: ' + df.get('fs_prefix', ''))
        click.echo('  Full path: ' + os.path.join(df.get('fs_prefix', ''), df['df']))
        click.echo('  Description: ' + df.get('description', ''))
        click.echo('')

@click.command(name='create')
@click.option('--df','-d', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--fs_prefix','-p', prompt='Filesystem prefix', help='Filesystem prefix' )
@click.pass_obj
def df_create(config, df, fs_prefix):
    """Create entry for dataset in database."""
    assnake_db_search = os.path.join(config['config']['assnake_db'], 'datasets/*')
    dfs = [d.split('/')[-1] for d in glob.glob(os.path.join(assnake_db_search))]
    print(df)
    print(fs_prefix)
    if df not in dfs:
        if os.path.isdir(os.path.join(fs_prefix, df)):
            df_info = {'df': df, 'fs_prefix': fs_prefix}
            os.makedirs(os.path.join(config['config']['assnake_db'], 'datasets/'+df), exist_ok=True)
            with open(os.path.join(config['config']['assnake_db'], 'datasets/'+df,'df_info.yaml'), 'w') as info_file:
                yaml.dump(df_info, info_file, default_flow_style=False)
            click.secho('Saved dataset ' + df + ' sucessfully!', fg='green')
    else:
        click.secho('Duplicate name!', fg='red')

@click.command(name ='info')
@click.option('--name','-d', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--preproc','-p', help='Show samples for preprocessing', required=False)
@click.pass_obj
def df_info(config, name, preproc):
    """View info for the specific dataset"""

    dfs = assnake.api.loaders.load_dfs_from_db('')
    df_info = dfs[name]
    click.echo(click.style(''*2 + df_info['df'] + ' '*2, fg='green', bold=True))
    click.echo('Filesystem prefix: ' + df_info.get('fs_prefix', ''))
    click.echo('Full path: ' + os.path.join(df_info.get('fs_prefix', ''), name))
    click.echo('Description: ' + df_info.get('description', ''))
    
    mg_samples_loc = os.path.join(config['config']['assnake_db'], 'datasets', df_info['df'], 'mg_samples.tsv')

    if os.path.isfile(mg_samples_loc):
        click.echo(tabulate(pd.read_csv(mg_samples_loc, sep='\t'), headers='keys', tablefmt='fancy_grid'))

    reads_dir = os.path.join(df_info['fs_prefix'], df_info['df'], 'reads/*')
    preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
    preprocs.sort()
    preprocessing = {}
    all_samples = []
    for p in preprocs:
        samples = assnake.api.sample_set.SampleSet(df_info['fs_prefix'], df_info['df'], p)
        # samples.add_samples(df_info['fs_prefix'], df_info['df'], p)
        all_samples += (list(samples.samples_pd['fs_name']))
        samples = samples.samples_pd[['fs_name', 'reads']].to_dict(orient='records')
        # click.secho(p + ' ' + str(len(samples)) + ' samples')
        preprocessing.update({p:samples})

    click.echo('\nTotal samples: ' + str(len(set(all_samples))) + '\n')


    for key, value in preprocessing.items():
        click.echo('Samples in ' + click.style(key, bold=True) + ': ' + str(len(value)))
    if preproc is not None:
        samples_pd = pd.DataFrame(preprocessing[preproc])
        # print(samples_pd)
        click.echo(tabulate(samples_pd.sort_values('reads'), headers='keys', tablefmt='fancy_grid'))
