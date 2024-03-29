import click, os, glob, json, zlib
import pandas as pd
from assnake.core.config import read_internal_config, read_assnake_instance_config, update_internal_config, fill_and_write_instance_config

from assnake.core.Result import Result
from assnake.core.snake_module import SnakeModule
import importlib

from assnake.utils.general import compute_crc32_of_dumped_dict
from pathlib import Path

from tabulate import tabulate

@click.command(name = 'list')
def show_installed_modules():
    """
    Show available assnake modules in your env
    """
    
    modules = SnakeModule.get_all_modules_as_dict()
    click.echo(modules)

@click.command(name = 'results')
def show_installed_results():
    """
    Show available results in your installation
    """
    
    results = Result.get_all_results_as_list()
    click.echo(results)

@click.command(name = 'redeploy')
@click.option('--module','-m', help='Snake Module to redeploy')
def redeploy_snake_module(module):
    """
    Redeploys Snake Module into database
    """
    sm = SnakeModule.get_all_modules_as_dict()[module]

    click.echo(sm)

@click.command(name = 'refresh-params')
def refresh_params():
    """
    Recomputes RCR32 hashes for parameters
    """
    tool = 'tmtic'

    instance_config = read_assnake_instance_config()
    if instance_config is not None:
        tool_params_dir = os.path.join(instance_config['assnake_db'], 'params/{tool}'.format(tool = tool))
        params_files = glob.glob(os.path.join(tool_params_dir, '*.json'))

        preset_table = []
        # Build preset table
        for preset_file_loc in params_files:
            splitted_name = os.path.basename(preset_file_loc).split('.')

            written_hash = 'NOT-PRESENT'
            if len(splitted_name) > 2:
                written_hash = os.path.basename(preset_file_loc).split('.')[1]
            
            computed_hash = compute_crc32_of_dumped_dict(preset_file_loc)
            ok = 'no'
            if computed_hash == written_hash:
                ok = 'yes'

            preset_table.append({
                'name': splitted_name[0],
                'ok': ok,
                'written_hash': written_hash,
                'computed_hash': computed_hash,
                'loc':  preset_file_loc,
            })
        preset_table = pd.DataFrame(preset_table)
        click.echo(tabulate(preset_table, headers='keys', tablefmt='fancy_grid'))

        not_ok = preset_table.loc[preset_table['ok'] == 'no']
        for preset in not_ok.to_dict(orient='records'):

            click.secho(' ==== ' + preset['name'] + ' ==== ', bold=True, fg='yellow')
            
            if preset['written_hash'] == 'NOT-PRESENT':
                warn_message = 'crc32 hash is not present in filename.'
                click.secho(warn_message, fg='yellow')

                # Check if there is another copy of this file
                is_unique = preset_table.loc[(preset_table['name'] == preset['name']) & (preset_table['computed_hash'] == preset['computed_hash'])]
                len_is_unique = len(is_unique)
                if len_is_unique > 1:
                    click.secho('We already have preset with this name and hash in database, so nothing to be done.', fg='green')
                elif len_is_unique == 1:
                    # Is it unique by name or by hash?
                    is_unique_name = preset_table.loc[(preset_table['name'] == preset['name'])]
                    is_unique_hash = is_unique.loc[(is_unique['computed_hash'] == is_unique['computed_hash'])]
                    if len(is_unique_name) == 1:
                        click.secho('This preset names is unique, just adding hash to the filename.', fg='green')
                    else:
                        click.secho('Name is not unique', fg='yellow')  
                        if len(is_unique_hash) == 1:
                            click.secho('It has duplicate name, but contents are different.\nThis might lead to confusion', fg='yellow')  
                            click.secho('It is STRONGLY advised to give your presets a prject specific unique name.', fg='yellow')  



def create_result_directory(result_name):
    """
    Creates a directory for the new result.
    """
    os.makedirs(result_name, exist_ok=True)
    return os.path.abspath(result_name)

def create_file_with_content(directory, filename, content):
    """
    Creates a file with the specified content.
    """
    with open(os.path.join(directory, filename), 'w') as file:
        file.write(content)

@click.command()
@click.option('--result-name', help='Name of your result. Will be used across system.')
@click.option('--description',  help='Brief description of the result')
@click.option('--result-type', default='preprocessing', help='Type of the result (e.g., preprocessing, analysis)')
@click.option('--input-type', default='illumina_sample', help='Input type for the result')
@click.option('--preset-file_format', default='yaml', help='Format of the preset file (yaml, json, etc.)')
def create_result(result_name, description, result_type, input_type, preset_file_format):
    """
    CLI command to create a new Result with a pre-populated template directory.
    """
    result_dir = create_result_directory(result_name)

    # Create result.py
    result_py_content = f'''
import os
from assnake.core.Result import Result

result = Result.from_location(name='{result_name}',
                              description='{description}',
                              result_type='{result_type}',
                              input_type='{input_type}',
                              with_presets=True,
                              preset_file_format='{preset_file_format}',
                              location=os.path.dirname(os.path.abspath(__file__)))
'''
    create_file_with_content(result_dir, 'result.py', result_py_content)

    # Create workflow.smk
    workflow_smk_content = '''# Add Snakemake workflow logic here'''
    create_file_with_content(result_dir, 'workflow.smk', workflow_smk_content)

    # Create wc_config.yaml
    wc_config_yaml_content = '''# Add wildcard configuration here'''
    create_file_with_content(result_dir, 'wc_config.yaml', wc_config_yaml_content)

    # Create wrapper.py
    wrapper_py_content = '''# Add wrapper logic here'''
    create_file_with_content(result_dir, 'wrapper.py', wrapper_py_content)

    click.echo(f'Result template created successfully at: {result_dir}')