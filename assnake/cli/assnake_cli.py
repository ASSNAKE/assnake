import sys, os, glob, yaml, shutil
import click
from assnake.core.Pipeline import Pipeline

import assnake.cli.commands.dataset_commands as dataset_commands # snakemake makes it slow
import assnake.cli.commands.config_commands as config_commands
import assnake.cli.commands.module_commands as module_commands
import assnake.cli.commands.pipeline_commands as pipeline_commands
from assnake.cli.commands.execute_commands import gather

from assnake.utils.general import read_yaml

from assnake.core.config import read_assnake_instance_config, read_internal_config, check_if_assnake_is_initialized

from pkg_resources import iter_entry_points 


#---------------------------------------------------------------------------------------
#                                    CLI  initialization
#---------------------------------------------------------------------------------------

@click.group()
@click.version_option()
@click.pass_context
def cli(ctx):
    """\b
   ___    ____   ____   _  __   ___    __ __   ____
  / _ |  / __/  / __/  / |/ /  / _ |  / //_/  / __/
 / __ | _\ \   _\ \   /    /  / __ | / ,<    / _/  
/_/ |_|/___/  /___/  /_/|_/  /_/ |_|/_/|_|  /___/


\b  
O---o Welcome to the Assnake, Traveller!
 O-o  If you found yorself here, 
  O   than you are exploring the land of the MICROBIOME.
 o-O  Here you can find weapons and spells
o---O that will help you to analyse your data.
O---o The tools that are presented here 
 O-o  are for metagenomics studies on ILLUMINA data.
  O   You can check quality and preprocess your data, 
 o-O  assemble it, map, bin the contigs, 
o---O make taxonomic annotations, functional annotations,
O---o and work with 16s rRNA data.
"""
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))

    # Load internal configuration but don't exit if instance configuration is missing
    internal_config = read_internal_config()
    instance_config = read_assnake_instance_config()
    if instance_config is not None:
        wc_config = read_yaml(os.path.join(dir_of_this_file, '../snake/wc_config.yaml'))
        ctx.obj = {'config': instance_config, 'wc_config': wc_config, 'requested_dfs': [], 'requests': [], 'sample_sets': [], 'requested_results': []}
    else:
        ctx.obj = {'config': None, 'wc_config': None, 'requested_dfs': [], 'requests': [], 'sample_sets': [], 'requested_results': []}

#---------------------------------------------------------------------------------------
#                                  assnake  INIT ***  group
#---------------------------------------------------------------------------------------

@cli.group(name='init')
def init_group():
    """Commands to initialize the ASSNAKE\n
    \bYou need to configure where assnake will store it's data and download databases.
    Assnake has internal configuration file 
    """
    pass


#---------------------------------------------------------------------------------------
#                                  assnake  DATASET ***  group
#---------------------------------------------------------------------------------------
@cli.group(chain=True)
def dataset():
    """Commands to work with datasets"""
    check_if_assnake_is_initialized()

dataset.add_command(dataset_commands.df_list)
dataset.add_command(dataset_commands.df_info)
dataset.add_command(dataset_commands.dataset_init)
dataset.add_command(dataset_commands.df_import_reads)


#---------------------------------------------------------------------------------------
#                                  assnake  RESULT ***  group
#---------------------------------------------------------------------------------------
from assnake.core.exceptions import InstanceConfigNotFound

@cli.group(chain=True, help='Used to request and run results')
def result():
    """Commands to analyze your data"""
    check_if_assnake_is_initialized()

try:
    for entry_point in iter_entry_points('assnake.plugins'):
        module_class = entry_point.load()
        module_name = entry_point.name  # Assuming this gives the module name

        for command in module_class.results:
            # Prefix command name with module name
            # TODO work on naming
            command_name = f"{module_name.replace('assnake-', '')}-{command.name}".replace('dada2-dada2-', 'dada2-').replace('core-preprocessing', 'qc')
            result.add_command(command.invocation_command, name=command_name)

    result.add_command(gather)

except InstanceConfigNotFound as e:
    click.secho(str(e), fg="red")




@cli.group(name = 'config')
def config_group():
    """Configuration related commands"""
    pass


config_group.add_command(config_commands.init_config)
config_group.add_command(config_commands.show_internal_config)
config_group.add_command(config_commands.show_instance_config)

@cli.group(name = 'module')
def module_group():
    """Commands to view and interact with assnake modules installed in current env"""
    pass

module_group.add_command(module_commands.show_installed_results)
module_group.add_command(module_commands.show_installed_modules)
module_group.add_command(module_commands.refresh_params)
module_group.add_command(module_commands.create_result)

@cli.group(name = 'pipeline')
def pipeline_group():
    """Commands to view and interact with assnake modules installed in current env"""
    pass

pipeline_group.add_command(pipeline_commands.run_pipeline)



def main():
    cli()

