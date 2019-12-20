import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import snakemake
from click.testing import CliRunner
from sys import argv

import assnake.commands.dataset_commands as dataset_commands

@click.group()
@click.version_option()
@click.pass_context
def cli(ctx):
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_loc = os.path.join(dir_of_this_file, '../snake/config.yml')

    if not os.path.isfile(config_loc):
        print("You need to init your installation! Iw won't take long. Just run assnake init start")
        exit()

    config = {}
    with open(config_loc, 'r') as stream:
        try:
            config = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)

    wc_config_loc = os.path.join(dir_of_this_file, '../snake/wc_config.yaml')
    wc_config = {}
    with open(wc_config_loc, 'r') as stream:
        try:
            wc_config = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)

    ctx.obj = {'config': config, 'wc_config': wc_config}
    pass 

@cli.group(name='init')
def init_group():
    """Commands to init the ASSNAKE"""
    pass

@click.command(name='start')
# @click.option('--config_location','-c', prompt='Desired location of the config file', help='Desired location of the config file' )
def init_start():
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_loc = os.path.join(dir_of_this_file, '../snake/config.yml')

    config = {}
    with open(config_loc, 'r') as stream:
        try:
            config = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)

    # if not os.path.isdir(config_location):
    #     print(config_location)
    #     os.makedirs(config_location, exist_ok=True)

    # with open(os.path.join(config_location, 'config.yml'), 'w') as config_file:
    #     written_config = yaml.dump(config, config_file)

    assnake_db = click.prompt('Please, enter desired location of the ASSNAKE database')
    config['assnake_db'] = assnake_db
    config['assnake_install_dir'] = os.path.join(os.path.dirname(dir_of_this_file), 'snakemake')
    
    with open(config_loc, 'w') as config_file:
        written_config = yaml.dump(config, config_file)

    # shutil.copyfile(config_loc, os.path.join(config_location, 'config.yml'))

init_group.add_command(init_start)

@cli.group(chain=True)
def dataset():
    """Commands to work with datasets"""
    pass

dataset.add_command(dataset_commands.df_list)
dataset.add_command(dataset_commands.df_info)
dataset.add_command(dataset_commands.df_create)

class ComplexCLI(click.MultiCommand):

    def list_commands(self, ctx):

        rv = []
        for filename in os.listdir('/data4/bio/fedorov/assnake_0.9.0/assnake/commands'):
            if filename.endswith('.py') and \
               filename.startswith('cmd_'):
                rv.append(filename[4:-3])


        for full in glob.glob('/data4/bio/fedorov/assnake_0.9.0/snake/modules/*/cmd_*.py'):
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