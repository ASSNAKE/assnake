import click, os, glob, json, zlib
import pandas as pd
from assnake.core.config import read_internal_config, read_assnake_instance_config, update_internal_config, fill_and_write_instance_config

from assnake.core.result import get_all_results_as_list
from assnake.core.snake_module import get_all_modules_as_dict
from assnake.utils.general import compute_crc32_of_dumped_dict
from pathlib import Path

from tabulate import tabulate

@click.command(name = 'list')
def show_installed_modules():
    """
    Show available assnake modules in your env
    """
    
    modules = get_all_modules_as_dict()
    click.echo(modules)

@click.command(name = 'results')
def show_installed_results():
    """
    Show available results in your installation
    """
    
    results = get_all_results_as_list()
    click.echo(results)

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

