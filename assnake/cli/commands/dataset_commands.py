#
# Perform assnake dataset [create/info/list] commands
#

import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
from assnake.api import fs_helpers
import importlib
from assnake.utils import pathizer, dict_norm_print
import snakemake



def show_av_dict(dfs):
    '''
    print out available datasets from dict of db's (assnake.api.loaders.load_dfs_from_db(''))
    '''
    avail_dfs = ''
    for i, item in enumerate(dfs.keys()):
        avail_dfs += '{}. {} \t'.format(i + 1, item)

    if avail_dfs == '':
        avail_dfs = 'NA'
    click.echo('Available: ' + avail_dfs)
    exit(2)


# ---------------------------------------------------------------------------------------
#                                     LIST
# ---------------------------------------------------------------------------------------
@click.command(name='list')
def df_list():
    """List datasets in database"""
    dfs = assnake.api.loaders.load_dfs_from_db('')
    if len(list(dfs.keys())) == 0:
        click.echo('No datasets in your system yet!\nYou can create one by running\n' +
                   click.style('  assnake dataset create  ', bg='blue', fg='white', bold=True))

    for df in dfs.values():
        df_name = df['df']
        click.echo(click.style('' * 2 + df_name + ' ' * 2, fg='green', bold=True))
        click.echo('  Filesystem prefix: ' + df.get('fs_prefix', ''))
        click.echo('  Full path: ' + os.path.join(df.get('fs_prefix', ''), df['df']))
        click.echo('  Description: ')
        dict_norm_print(df.get('description', ''))
        click.echo('')


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# ---------------------------------------------------------------------------------------
#                                   CREATE
# ---------------------------------------------------------------------------------------
@click.command(name='create')
@click.option('--df', '-d', prompt='Name of the dataset', help='Name of the dataset')
@click.option('--fs_prefix', '-f', prompt='Filesystem prefix', help='Filesystem prefix. This MUST be ABSOLUTE path')
@click.pass_obj
def df_create(config, df, fs_prefix):
    there_is_prpties = False
    """Register your dataset inside ASSNAKE!\n
    You can use it in interactive mode."""
    if fs_prefix[0] != '/':
        click.echo('THIS IS NOT AN ABSOLUTE PATH!!!')
        exit()
    assnake_db_search = os.path.join(config['config']['assnake_db'], 'datasets/*')
    dfs = [d.split('/')[-1] for d in glob.glob(os.path.join(assnake_db_search))]

    prpts = dict()
    while click.confirm('Add property', abort=False):
        if not there_is_prpties:
            there_is_prpties = True
        prpts.update({click.prompt('Name of property:'): click.prompt('Value of property:')})

    if df not in dfs:
        if os.path.isdir(os.path.join(fs_prefix, df)):
            df_info = {'df': df, 'fs_prefix': fs_prefix, 'description': prpts}
            os.makedirs(os.path.join(config['config']['assnake_db'], 'datasets/' + df), exist_ok=True)
            with open(os.path.join(config['config']['assnake_db'], 'datasets/' + df, 'df_info.yaml'), 'w') as info_file:
                yaml.dump(df_info, info_file, default_flow_style=False)
            click.secho('Saved dataset ' + df + ' sucessfully!', fg='green')
        else:
            click.secho('We were unable to find directory on the path')
    else:
        click.secho('Duplicate name!', fg='red')
        if there_is_prpties:
            click.echo('Description updated')
            with open(os.path.join(config['config']['assnake_db'], 'datasets/' + df, 'df_info.yaml'),
                      'r') as info_file_old:
                df_info_old = yaml.load(info_file_old)
            try:
                df_info_old['description'].update(prpts)
            except KeyError as e:
                df_info_old['description'] = prpts
            with open(os.path.join(config['config']['assnake_db'], 'datasets/' + df, 'df_info.yaml'), 'w') as info_file:
                yaml.dump(df_info_old, info_file, default_flow_style=False)


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# ---------------------------------------------------------------------------------------
#                                   INFO
# ---------------------------------------------------------------------------------------
@click.command(name='info')
@click.option('--name', '-d', prompt='Name of the dataset', help='Name of the dataset')
@click.option('--preproc', '-p', help='Show samples for preprocessing', required=False)
@click.pass_obj
def df_info(config, name, preproc):
    """View info for the specific dataset"""
    dfs = assnake.api.loaders.load_dfs_from_db('')
    try:
        df_info = dfs[name]
    except Exception as e:
        click.echo('Can`t reach database with such name')
        show_av_dict(dfs)
    click.echo(click.style('' * 2 + df_info['df'] + ' ' * 2, fg='green', bold=True))
    click.echo('Filesystem prefix: ' + df_info.get('fs_prefix', ''))
    click.echo('Full path: ' + os.path.join(df_info.get('fs_prefix', ''), name))
    click.echo('Description: ')
    dict_norm_print(df_info.get('description', ''))

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
        preprocessing.update({p: samples})

    click.echo('\nTotal samples: ' + str(len(set(all_samples))) + '\n')

    for key, value in preprocessing.items():
        click.echo('Samples in ' + click.style(key, bold=True) + ': ' + str(len(value)))
    if preproc is not None:
        samples_pd = pd.DataFrame(preprocessing[preproc])
        # print(samples_pd)
        click.echo(tabulate(samples_pd.sort_values('reads'), headers='keys', tablefmt='fancy_grid'))


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# ---------------------------------------------------------------------------------------
#                                   IMPORT-READS
# ---------------------------------------------------------------------------------------
# DONE decide if we need either d and t or proceed both arguments as one and automatically choose path or not
@click.command(name='import-reads')
@click.option('--source', '-r', prompt='Location of the source dataset.',
              help='Location of the source dataset. Absolute path required.', type=click.Path())
@click.option('--dataset', '-d', help='Assnake dataset name', required=False)
@click.option('--target', '-t', help='Location of the target directory. Absolute path required.', required=False,
              type=click.Path())
@click.option('--sample_set', '-s', help='Comma-divided list of samples of interest', required=False)
@click.option('--sample_list', '-l', help='Absolute path to file with line by line samples of interest', required=False,
              type=click.Path())
@click.option('--copy', help='If is set, hard copying will be used instead of symbolic links ', is_flag=True)
@click.pass_obj
def df_import_reads(config, source, dataset, target, sample_set, sample_list, copy):
    """
    Import reads from collective file to individual files directory
    """
    arg_d = not bool(dataset is None)
    arg_t = not bool(target is None)
    arg_s = not bool(sample_set is None)
    arg_l = not bool(sample_list is None)


    if arg_l & arg_s:
        click.secho('Collision tends to be observed. Please, specify either list of samples in prompt or in file',
                    err=True)
        exit(1)
    if not (arg_d ^ arg_t):
        click.secho('Please, specify either database (-d) or absolute path (-t)', err=True)
        exit(1)
    importlib.reload(fs_helpers)

    if arg_d:
        try:
            df_info = assnake.api.loaders.load_df_from_db(dataset)
        except Exception as e:
            dfs = assnake.api.loaders.load_dfs_from_db('')
            click.echo('Can`t reach database with such name', err=True)
            show_av_dict(dfs)
        target = '{}/{}/reads/raw'.format(df_info['fs_prefix'], df_info['df'])
    else:
        target = pathizer(target)
        if not os.path.exists(target):
            click.secho("Provided sample-list file couldn't be detected", err=True)
            exit(2)
    print(target)

    # def rename(sample):
    #     new_name_wc = 'Rst{}_{}_B2'
    #     splitted = sample.split('_')
    #     splitted[2] = (1 if len(splitted) == 3 else a)
    #     return new_name_wc.format(splitted[1], splitted[2])

    dicts = fs_helpers.get_samples_from_dir(source)
    sample_names = {d['sample_name'] for d in dicts}
    if arg_s:
        samples_of_interest = sample_set.split(',')
    elif arg_l:
        sample_list = pathizer(sample_list)
        if os.path.exists(sample_list):
            with open(sample_list, 'r') as ls_file:
                samples_of_interest = ls_file.readlines()
        else:
            click.secho("Provided sample-list file couldn't be detected", err=True)

    else:
        samples_of_interest = sample_names
    if (arg_s or arg_l) and (len(set(samples_of_interest) - {d['sample_name'] for d in dicts}) != 0):
        click.secho("Warning! Samples, been specified, are not in reads file", err=True)
        click.confirm('Continue?', abort=True)
    if copy:
        click.secho("Please, keep in mind, that hard copying may take plenty of time")
    for d in dicts:
        if d['sample_name'] in samples_of_interest:
            fs_helpers.create_links(source, target, d, hard=copy)