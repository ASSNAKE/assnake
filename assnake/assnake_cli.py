import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import snakemake
from click.testing import CliRunner
from sys import argv

import assnake.commands.dataset_commands as dataset_commands
import assnake.commands.init as commands_init
import configparser
import shutil
import pprint

from assnake.utils import read_yaml

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
    config_internal = configparser.ConfigParser()
    config_internal.read(os.path.join(dir_of_this_file, './config_internal.ini'))
    config_loc = config_internal['GENERAL']['config_loc']

    if not os.path.isfile(config_loc):
        print("You need to init your installation! Iw won't take long. Just run assnake init start")
        exit()
    else:
        # click.echo('We found config at ' + config_loc)
        pass

    config = read_yaml(config_loc)
    wc_config_loc = os.path.join(dir_of_this_file, '../snake/wc_config.yaml')
    wc_config = read_yaml(wc_config_loc)

    ctx.obj = {'config': config, 'wc_config': wc_config}
    pass 



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
    pass

dataset.add_command(dataset_commands.df_list)
dataset.add_command(dataset_commands.df_info)
dataset.add_command(dataset_commands.df_create)

class ComplexCLI(click.MultiCommand):

    def list_commands(self, ctx):
        dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
        dir_for_commands = os.path.join(dir_of_this_file, 'commands')
        dir_for_snake_modules = os.path.join(dir_of_this_file, '../snake/modules/*/cmd_*.py')

        rv = []
        for filename in os.listdir(dir_for_commands):
            if filename.endswith('.py') and \
               filename.startswith('cmd_'):
                rv.append(filename[4:-3])


        for full in glob.glob(dir_for_snake_modules):
            filename = full.split('/')[-1]
            module = full.split('/')[-2]
            if filename.endswith('.py') and \
               filename.startswith('cmd_'):
                rv.append(module +'/'+filename[4:-3])

        print( rv)
        # rv.sort()
        return rv

    def get_command(self, ctx, name):
        if '/' not in name and name in ['construct_sample_set_for_assembly', 'run', 'status']:
            try:
                mod = __import__('assnake.commands.cmd_' + name,
                                None, None, ['cli'])
            except ImportError as error:
                print('import_error')
                print(error.__class__.__name__ + ": ", name)
                return None
            return mod.cli
        else:
            try:
                module = name.split('/')[0]
                command_name = name.split('/')[-1]
                mod_wc = 'snake.modules.{module}.cmd_{command_name}'
                mod = __import__(mod_wc.format(module = module, command_name = command_name),
                                None, None, ['cli'])
            except ImportError as error:
                print('=import_error=')
                print(error.__class__.__name__ + ": ", name)

                return None
            return mod.cli
        

@cli.group()
def result():
    pass

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

@click.group(cls=ComplexCLI, chain = True, help = 'Used to request and run results')
@click.pass_obj
def request(config):
    pass

result.add_command(request)



def main():
    # TODO Here we should manually parse and check that if we request result `run` command is last
    # print('======IN MAIN===========')
    # print(argv[:])
    # try:
    #     print(argv[:].index('run'))
    # except:
    #     print('no run in index')
    cli()


    #     1 - Quality Control
    #     FastQC - Check quality of your reads
    #     Trimmomatic 
    #     Remove human reads with bbmap
    #     Filter reads by minimum length
    #     MultiQC
    # 2 - Assembly
    #     Megahit
    #     Spades (need to include parameters parser)
    # 3 - Taxonomic Annotation
    #     Metaphlan2
    #     Centrifuge
    #     Kraken (Rewrite needed)
    # 4 - Functional Annotation
    #     Humann2 (Conflict with current metaphlan version!!)
    # 5 - Binning
    #     MaxBin 2
    #     MetaBat 2
    # 6 - 16S
    #     Dada2
    #     Picrust 2

# ╔═╗╔═╗╔═╗╔╗╔╔═╗╦╔═╔═╗
# ╠═╣╚═╗╚═╗║║║╠═╣╠╩╗║╣ 
# ╩ ╩╚═╝╚═╝╝╚╝╩ ╩╩ ╩╚═╝
# \b