import click, sys, os, glob, yaml, shutil
import pandas as pd
from tabulate import tabulate
from assnake.utils import fs_helpers
from assnake.utils.general import pathizer, dict_norm_print, download_from_url
from zipfile import ZipFile
from pathlib import Path

from assnake.core.Dataset import Dataset
from assnake.core.SampleContainerSet import SampleContainerSet
from yaml import dump

import os
from pathlib import Path
from shutil import copy2

# ---------------------------------------------------------------------------------------
#                                     LIST
# ---------------------------------------------------------------------------------------
@click.command(name='list')
def df_list():
    """List datasets in database"""
    dfs = Dataset.list_in_db()
    if len(list(dfs.keys())) == 0:
        click.echo('No datasets in your system yet!\nYou can create one by running\n' +
                   click.style('  assnake dataset create  ', bg='blue', fg='white', bold=True))

    for df in dfs.values():
        df_name = df['df']
        click.echo(click.style('' * 2 + df_name + ' ' * 2, fg='green', bold=True))
        # click.echo('  Filesystem prefix: ' + df.get('fs_prefix', ''))
        click.echo('  Full path: ' + os.path.join(df.get('fs_prefix', ''), df['df']))
        # click.echo('  Description: ')
        # dict_norm_print(df.get('description', ''))
        click.echo('')





def check_dataset_availability(dataset_name, datasets_in_db):
    """
    Check if the dataset name already exists in the database.
    """
    if dataset_name in datasets_in_db:
        click.secho('Dataset with the name already exists. Aborting.', fg='red')
        raise click.Abort()

def validate_dataset_path(fs_prefix, dataset_name, allow_non_empty=False):
    """
    Validate and construct the full path for the dataset.
    """
    full_dataset_path = os.path.join(fs_prefix, dataset_name)
    if os.path.isdir(full_dataset_path):
        if not os.listdir(full_dataset_path) or allow_non_empty:
            return full_dataset_path
        else:
            click.secho('Dataset path exists and is not empty. Aborting.', fg='red')
            raise click.Abort()
    else:
        os.makedirs(full_dataset_path, exist_ok=True)
        return full_dataset_path

def create_symlink_to_dataset(assnake_db, dataset_name, full_dataset_path):
    """
    Create a symbolic link in the Assnake database pointing to the dataset.
    """
    symlink_path = os.path.join(assnake_db, 'datasets', dataset_name)
    if not os.path.islink(symlink_path):
        os.symlink(full_dataset_path, symlink_path, target_is_directory=True)

def create_dataset_info_file(dataset_path, dataset_info):
    """
    Create and write the dataset information file.
    """
    with open(os.path.join(dataset_path, 'df_info.yaml'), 'w') as info_file:
        dump(dataset_info, info_file, default_flow_style=False)

# ---------------------------------------------------------------------------------------
#                                   INIT
# ---------------------------------------------------------------------------------------

@click.command(name='init')
@click.option('--data-type', '-t', type=click.Choice(['METAGENOMIC_16s', 'METAGENOMIC_WGS', 'VIROME', 'RNA_SEQ'], case_sensitive=False))
@click.option('--dataset-name', '-d', required=False)
@click.option('--first-preprocessing-name', default='raw')
@click.pass_obj
def dataset_init(config, data_type, dataset_name, first_preprocessing_name):
    """
    Initialize a dataset within the Assnake framework.
    """
    cwd = os.getcwd()
    datasets_in_db = Dataset.list_in_db()

    if not dataset_name:
        dataset_name = os.path.basename(cwd)
    check_dataset_availability(dataset_name, datasets_in_db)

    fs_prefix = cwd if not dataset_name else os.path.dirname(cwd)
    full_dataset_path = validate_dataset_path(fs_prefix, dataset_name)
    
    os.makedirs(os.path.join(full_dataset_path, 'reads/raw'))

    dataset_info = {'df': dataset_name, 'fs_prefix': fs_prefix, 'description': {}, 'data_type': data_type}
    create_symlink_to_dataset(config['config']['assnake_db'], dataset_name, full_dataset_path)
    create_dataset_info_file(full_dataset_path, dataset_info)

    click.secho(f"Dataset '{dataset_name}' initialized successfully.", fg='green')

# ---------------------------------------------------------------------------------------
#                                   INFO
# ---------------------------------------------------------------------------------------
@click.command(name='info')
@click.option('--df', '-d', help='Name of the dataset', required=False)
@click.option('--preproc', '-p', help='Show samples for preprocessing', required=False)
@click.argument('df_arg', required=False)
@click.pass_obj
def df_info(config, df, preproc, df_arg):
    """View info for the specific dataset
        Usage: assnake dataset info [dataset] or -d [dataset] ...

    """
    # df argument/option logic
    if not (bool(df is None) ^ bool(df_arg is None)):
        click.echo('Please, specify dataset either as option or argument')
        df = click.prompt('Type the name in')
    if df is None:
        df = df_arg
    
    dataset = Dataset(df)
    # Formatting a dataset name as a header
    header = f'==  {dataset.dataset_name}  =='
    styled_header = click.style(header, fg='green', bold=True)
    click.echo(styled_header)

    # Displaying dataset information
    click.echo(str(dataset))
    # Now create a SampleContainerSet with specific parameters
    sample_container_set = SampleContainerSet(dataset, preproc='raw')

    for sample_container in sample_container_set.sample_containers:
        click.echo(sample_container)

    return



# ---------------------------------------------------------------------------------------
#                                   IMPORT-READS
# ---------------------------------------------------------------------------------------




def get_target_directory(dataset_name, target, config):
    """
    Determine the target directory for importing reads.
    """
    if dataset_name:
        try:
            dataset_info = Dataset(dataset_name)
            return os.path.join(dataset_info.full_path, 'reads', 'raw')
        except :
            available_datasets = Dataset.list_in_db()
            click.secho(f"Cannot find dataset: {dataset_name}", fg='red')
            # Show available datasets
            for i, ds in enumerate(available_datasets.keys(), 1):
                click.echo(f'{i}. {ds}')
            raise click.Abort()
    else:
        target = str(Path(target).resolve())
        if not os.path.exists(target):
            click.secho("Target directory does not exist: " + target, fg='red')
            raise click.Abort()
        return target

def modify_sample_name(name, rename_method):
    """
    Modify the sample name based on the chosen rename method.
    """
    if rename_method == 'removeSending':
        return '_'.join(name.replace('-', '_').split('_')[:-1])
    else:  # Default to 'replace-'
        return name.replace('-', '_')

@click.command(name='import-reads')
@click.option('--reads-dir', '-r', prompt='Location of folder with read files', type=click.Path(exists=True), help='Directory containing read files')
@click.option('--dataset', '-d', help='Name of the Assnake dataset', required=False)
@click.option('--rename-method', help='Method to rename samples', type=click.Choice(['replace-', 'removeSending'], case_sensitive=False), default='replace-')
@click.option('--target', '-t', help='Target directory for imported reads', type=click.Path(), required=False)
@click.option('--copy', is_flag=True, help='Use hard copying instead of symbolic links')
@click.pass_obj
def df_import_reads(config, reads_dir, dataset, rename_method, target, copy):
    """
    Import reads into an Assnake dataset.
    """
    # Validate and get target directory
    target_dir = get_target_directory(dataset, target, config) if dataset or target else None
    if not target_dir:
        click.secho("Target directory is not specified. Aborting.", fg='red')
        raise click.Abort()

    # Ensure target directory exists
    os.makedirs(target_dir, exist_ok=True)

    # Modify sample names and import reads
    modify_name = lambda name: modify_sample_name(name, rename_method)
    samples_in_run = fs_helpers.get_samples_from_dir(reads_dir, modify_name)
    if samples_in_run.empty:
        click.secho('No reads found in the specified directory.', fg='yellow')
        return

    samples_in_run['df_sample'] = samples_in_run['modified_name']
    fs_helpers.create_links(target_dir, samples_in_run, hard=copy)
    click.secho("Reads imported successfully.", fg='green')