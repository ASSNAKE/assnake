import click, os

from assnake.core.config import read_internal_config, read_assnake_instance_config, update_internal_config, fill_and_write_instance_config

from assnake.core.result import get_all_results_as_list
from assnake.core.snake_module import get_all_modules_as_dict

from pathlib import Path

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