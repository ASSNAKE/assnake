import sys, os, glob, yaml, shutil
import click

import assnake.cli.commands.dataset_commands as dataset_commands # snakemake makes it slow
import assnake.cli.commands.init_commands as init_commands
from assnake.cli.commands.execute_commands import gather

import assnake.utils
from assnake.utils import read_assnake_instance_config, read_internal_config, read_yaml, check_if_assnake_is_initialized
from assnake.cli.cli_utils import sample_set_construction_options, add_options
from pkg_resources import iter_entry_points 
from assnake.core.sample_set import generic_command_individual_samples, generate_result_list


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
O---o
 O-o  and work with 16s rRNA data.
  O
 o-O
o---O

\b
Please, feel free to check out all the commands. 
Please, read the help messages and if you have any questions, 
write me directly on my email fedorov.de@gmail.com, 
or create an issue or PR on GitHub.
\b
Start by initializing ASSNAKE with
assnake init command
\b"""


    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))

    instance_config = read_assnake_instance_config()

    if instance_config is not None:
        wc_config = read_yaml(os.path.join(dir_of_this_file, '../snake/wc_config.yaml'))
        ctx.obj = {'config': instance_config, 'wc_config': wc_config, 'requested_dfs': [], 'requests': [], 'sample_sets': [], 'requested_results': []}


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

init_group.add_command(init_commands.init_start)

#---------------------------------------------------------------------------------------
#                                  assnake  DATASET ***  group
#---------------------------------------------------------------------------------------
@cli.group(chain=True)
def dataset():
    """Commands to work with datasets"""
    check_if_assnake_is_initialized()

dataset.add_command(dataset_commands.df_list)
dataset.add_command(dataset_commands.df_info)
dataset.add_command(dataset_commands.df_init)
dataset.add_command(dataset_commands.df_create)
dataset.add_command(dataset_commands.df_import_reads)
dataset.add_command(dataset_commands.df_delete)
dataset.add_command(dataset_commands.rescan_dataset)

#---------------------------------------------------------------------------------------
#                                  assnake  RESULT ***  group
#---------------------------------------------------------------------------------------
@cli.group(chain = True, help = 'Used to request and run results')
def result():
    """Commands to analyze your data"""
    check_if_assnake_is_initialized()


for entry_point in iter_entry_points('assnake.plugins'):
    module_class = entry_point.load()
    for cmd in module_class.invocation_commands:
        result.add_command(cmd)
    for cmd in module_class.initialization_commands:
        init_group.add_command(cmd)
    for snakeresult in module_class.results:
        result.add_command(snakeresult.invocation_command)

result.add_command(gather)


@click.command('sample-set', short_help='Filter and trim your reads with dada2 trimmer')
@add_options(sample_set_construction_options)
@click.option('--params', help='Parameters to use', default='def', type=click.STRING )
@click.option('--message', '-m', multiple=True)

@click.pass_obj
def request_sample_set(config, message, **kwargs):
    click.echo('\n'.join(message))
    sample_set, sample_set_name = generic_command_individual_samples(config,  **kwargs)
    config['sample_sets'].append(sample_set.samples_pd.copy())
result.add_command(request_sample_set)


@cli.group(name = 'config')
def config_group():
    """Configuration related commands"""
    pass

@config_group.command(name = 'show-internal')
def show_internal_config():
    """
    Show your current internal configuration
    """
    # click.echo('CURRENT CONFIG LOCATION')
    print(assnake.utils.read_internal_config())

@config_group.command(name = 'show-instance')
def show_internal_config():
    """
    Show your current instance configuration
    """
    # click.echo('CURRENT CONFIG LOCATION')
    print(assnake.utils.read_assnake_instance_config())

def main():
    cli()

