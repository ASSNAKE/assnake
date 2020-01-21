import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import snakemake
from click.testing import CliRunner
from sys import argv

import assnake.cli.commands.dataset_commands as dataset_commands
import assnake.cli.commands.init as commands_init
from assnake.cli.commands.cmd_run import run
import configparser
import shutil
import pprint

from assnake.utils import read_yaml

import pkg_resources



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
\b
Here is it how it goes.
Somewhere on your filesystem you create a folder, and put your reads inside the ./<your_folder>/reads/raw folder.
<your_folder> is also the name of the Dataset, so choose wisely!
Than register yor dataset in the assnake with
assnake dataset create

    
    """
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_loc = assnake.utils.get_config_loc()

    if not os.path.isfile(config_loc):
        pass
    else:
        config = read_yaml(config_loc)
        wc_config = read_yaml(os.path.join(dir_of_this_file, '../snake/wc_config.yaml'))

        ctx.obj = {'config': config, 'wc_config': wc_config}



@cli.group(name='init')
def init_group():
    """Commands to initialize the ASSNAKE\n
    \bYou need to configure where assnake will store it's data and download databases.
    Assnake has internal configuration file 
    """
    pass
init_group.add_command(commands_init.init_start)

@cli.group(chain=True)
def dataset():
    """Commands to work with datasets"""
    assnake.utils.check_if_assnake_is_initialized()

dataset.add_command(dataset_commands.df_list)
dataset.add_command(dataset_commands.df_info)
dataset.add_command(dataset_commands.df_create)

@cli.group()
def result():
    assnake.utils.check_if_assnake_is_initialized()

@result.command('prepare_set_for_assembly')
@click.option('--df', '-d')
@click.option('--preproc','-p', help='Preprocessing to use' )
@click.option('--samples-to-add','-s', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING )
@click.option('--set_name','-n', help='Name of your sample set' )


def prepare_set_for_assembly(df, preproc, samples_to_add,set_name):
    print(df)
    print(preproc)
    print(samples_to_add)
    samples_to_add = [] if samples_to_add == '' else [c.strip() for c in samples_to_add.split(',')]
    df = assnake.api.loaders.load_df_from_db(df)
    ss = assnake.api.sample_set.SampleSet(df['fs_prefix'], df['df'], preproc, samples_to_add=samples_to_add)

    click.echo(tabulate(ss.samples_pd[['fs_name', 'reads', 'preproc', 'df']].sort_values('reads'), 
        headers='keys', tablefmt='fancy_grid'))

    set_dir = os.path.join(df['fs_prefix'], df['df'], 'assembly', set_name)
    os.makedirs(set_dir, exist_ok=True)

    set_loc = os.path.join(set_dir, 'sample_set.tsv')
    ss.samples_pd[['df', 'preproc', 'fs_name']].to_csv(set_loc, sep='\t', index=False)

@click.group(chain = True, help = 'Used to request and run results')
@click.pass_obj
def request(config):
    pass

result.add_command(request)

discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point in pkg_resources.iter_entry_points('assnake.plugins')
}
for module_name, module_class in discovered_plugins.items():
    for cmd in module_class.invocation_commands:
        request.add_command(cmd)
request.add_command(run)


def main():
    # TODO Here we should manually parse and check that if we request result `run` command is last
    # print('======IN MAIN===========')
    # print(argv[:])
    # try:
    #     print(argv[:].index('run'))
    # except:
    #     print('no run in index')
    cli()

