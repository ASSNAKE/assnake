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

    # print(config_loc)
    if not os.path.isfile(config_loc):
        print("You need to init your installation! Iw won't take long. Just run assnake init start")
        exit()
    # else:
    #     print("we found default config file")

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
    pass #Entry Point

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

@cli.group()
def result():
    """Commands to work with results"""
    pass

example_res = ['fastqc', 'count', 'metaphlan2']

@click.command()
@click.option('--df','-d', help='Name of the dataset' )
@click.option('--preproc','-p', help='Preprocessing to use' )
@click.option('--samples-to-add','-r', 
                help='Samples from dataset to process', 
                default='', 
                metavar='<samples_to_add>', 
                type=click.STRING )

@click.option('--results','-r', 
                prompt='Results you want', 
                help='Results you want', 
                default=','.join(example_res), 
                show_default=True,
                metavar='<results>', 
                type=click.STRING )

@click.option('--params', help='Parameters id to use', required=False, default = 'def')
@click.option('--list-name', help='Name of sample set if required by rule')

@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--run/--no-run', default=False)
@click.pass_obj
def request(config, df, preproc, samples_to_add, results, params, list_name,   threads, jobs, run):
    """Request result for your samples.\nYou need to provide info defining sample set, or run it from {prefix}/{df}/reads/{preproc} directory"""
    if samples_to_add == '':
        samples_to_add = []
    else:
        samples_to_add = [c.strip() for c in samples_to_add.split(',')]
    # print(samples_to_add)
    ss = None
    if df is not None:
        df = assnake.api.loaders.load_df_from_db(df)
        # assnake.api.sample_set.SampleSet(df['fs_prefix'], df['df'], preproc, samples_to_add=samples_to_add)
        ss = assnake.api.sample_set.SampleSet(df['fs_prefix'], df['df'], preproc, samples_to_add=samples_to_add)
        click.echo(tabulate(ss.samples_pd[['fs_name', 'reads', 'preproc']].sort_values('reads'), 
                headers='keys', tablefmt='fancy_grid'))

    res_list = []
    # split results by ',' and remove whitespace
    results = [c.strip() for c in results.split(',')]
    for result in results:
        if result == 'multiqc':
            ss.prepare_fastqc_list_multiqc('R1', preproc)
            ss.prepare_fastqc_list_multiqc('R2', preproc)
            multiqc_report_wc = config['wc_config']['multiqc_report']

            res_list = [
                multiqc_report_wc.format(
                    prefix = df['fs_prefix'],
                    df = df['df'],
                    sample_set = preproc,
                    strand = strand,
                )
                for strand in ['R1', 'R2']
            ]
        elif result == 'megahit':
            # check for prepared sample sets
            prepared_sets = glob.glob(os.path.join(df['fs_prefix'],df['df'],'assembly/*/*/sample_set.tsv'))
            click.echo('We found ' + str(len(prepared_sets)) + ' prepared sets:')
            for i, pset in enumerate(prepared_sets):
                name = pset.split('/')[-2]
                run_info = pset.split('/')[-3]
                click.echo(click.style(str(i), bold = True) + ' ' + run_info + ' ' + name)

            selected_sets = []
            sel_sets_str = click.prompt('Please, enter sets you want to assemble, comma separated. You can also type all, to select all sets', type=str)
            sel_sets_str = sel_sets_str.replace(' ', '')
            if sel_sets_str != 'all':
                selected_sets = [int(s) for s in sel_sets_str.split(',')] 
            click.echo(click.style('Selected sets: ' + str(set(selected_sets)), fg='green'))
            min_len = click.prompt('Please, enter minimum contig lenth', type=int)

            for ss in selected_sets:
                res_list += [prepared_sets[ss].replace('sample_set.tsv', 'final_contigs__{min_len}.fa'.format(min_len=min_len))]
        elif result=='metabat2':
            mb2 = '{fs_prefix}/{df}/metabat2/bwa__0.7.17__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/metabat2.done'

            prepared_sets = glob.glob(os.path.join(df['fs_prefix'],df['df'],'assembly/*/*/sample_set.tsv'))
            click.echo('We found ' + str(len(prepared_sets)) + ' prepared sets:')
            for i, pset in enumerate(prepared_sets):
                name = pset.split('/')[-2]
                run_info = pset.split('/')[-3]
                click.echo(click.style(str(i), bold = True) + ' ' + run_info + ' ' + name)

            selected_sets = []
            sel_sets_str = click.prompt('Please, enter sets you want to assemble, comma separated. You can also type all, to select all sets', type=str)
            sel_sets_str = sel_sets_str.replace(' ', '')
            if sel_sets_str != 'all':
                selected_sets = [int(s) for s in sel_sets_str.split(',')] 
            else: 
                selected_sets = list(range(0,len(prepared_sets)))
            click.echo(click.style('Selected sets: ' + str(set(selected_sets)), fg='green'))
            min_len = click.prompt('Please, enter minimum contig lenth', type=int)

            for ss in selected_sets:
                name = prepared_sets[ss].split('/')[-2]
                res_list += [mb2.format(fs_prefix=df['fs_prefix'],df=df['df'], params='def',sample_set=name, mod=min_len)]
        elif result=='maxbin2':
            mb2 = '{fs_prefix}/{df}/maxbin2/bwa__0.7.17__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/maxbin2.done'

            prepared_sets = glob.glob(os.path.join(df['fs_prefix'],df['df'],'assembly/*/*/sample_set.tsv'))
            click.echo('We found ' + str(len(prepared_sets)) + ' prepared sets:')
            for i, pset in enumerate(prepared_sets):
                name = pset.split('/')[-2]
                run_info = pset.split('/')[-3]
                click.echo(click.style(str(i), bold = True) + ' ' + run_info + ' ' + name)

            selected_sets = []
            sel_sets_str = click.prompt('Please, enter sets you want to assemble, comma separated. You can also type all, to select all sets', type=str)
            sel_sets_str = sel_sets_str.replace(' ', '')
            if sel_sets_str != 'all':
                selected_sets = [int(s) for s in sel_sets_str.split(',')] 
            else: 
                selected_sets = list(range(0,len(prepared_sets)))
            click.echo(click.style('Selected sets: ' + str(set(selected_sets)), fg='green'))
            min_len = click.prompt('Please, enter minimum contig lenth', type=int)

            for ss in selected_sets:
                name = prepared_sets[ss].split('/')[-2]
                res_list += [mb2.format(fs_prefix=df['fs_prefix'],df=df['df'], params='def',sample_set=name, mod=min_len)]
        elif result=='mp2':
            res_list = ss.get_locs_for_result(result)
        elif result=='count':
            res_list += ss.get_locs_for_result(result)
        elif result=='fastqc':
            res_list += ss.get_locs_for_result(result)
        elif result=='trimmomatic':
            res_list += ss.get_locs_for_result(result, params=params)
        elif result == 'dada2-filter-and-trim':
            res_list = ss.get_locs_for_result(result, params=params)
        elif result=='dada2-learn-errors':
            # ss.prepare_dada2_sample_list()
            if list_name is None:
                list_name = click.prompt('Please enter name for your dada2 error profile')
                print(list_name)
            if df is not None:
                ss.prepare_dada2_sample_list(list_name)
                res_list= [(config['wc_config']['dada2_err_wc'].format(
                        dada2_dir = config['config']['dada2_dir'],
                        sample_set = list_name,
                        strand = strand,
                        params=params
                    )) for strand in ['R1', 'R2']]
            else:
                res_list= [(config['wc_config']['dada2_err_wc'].format(
                        dada2_dir = config['config']['dada2_dir'],
                        sample_set = list_name,
                        strand = strand,
                        params=params
                    )) for strand in ['R1', 'R2']]
        

    
    curr_dir = os.path.abspath(os.path.dirname(__file__))
    click.echo('Current dir: ' + curr_dir)
    click.echo('Number of jobs torun in parallel: ' + str(jobs))

    status = snakemake.snakemake(os.path.join(curr_dir, './bin/snake/base.py'), 
        # config = config['wc_config'],    
        targets=res_list, 
        printshellcmds=True,
        dryrun=not run, 
        configfiles=[os.path.join(curr_dir, '../snake/config.yml')],
        drmaa_log_dir = config['config']['drmaa_log_dir'],
        use_conda = True,
        latency_wait = 120,
        conda_prefix = config['config']['conda_dir'],
        drmaa=' -V -S /bin/bash -pe make {threads}'.format(threads=threads),
        cores=jobs, nodes=jobs
        )

@click.command(name='list')
@click.pass_obj
def res_list(config):
    click.echo('results')

result.add_command(request)
result.add_command(res_list)


class ComplexCLI(click.MultiCommand):

    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir('/data4/bio/fedorov/assnake_0.9.0/assnake/commands'):
            if filename.endswith('.py') and \
               filename.startswith('cmd_'):
                rv.append(filename[4:-3])
        for filename in os.listdir('/data4/bio/fedorov/assnake_0.9.0/snake/modules/fastqc'):
            if filename.endswith('.py') and \
               filename.startswith('cmd_'):
                rv.append(filename[4:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        try:
            if sys.version_info[0] == 2:
                name = name.encode('ascii', 'replace')
            mod = __import__('assnake.commands.cmd_' + name,
                             None, None, ['cli'])
        except ImportError as error:
            print('import_error')
            print(error.__class__.__name__ + ": ")
            print(name)
            mod = __import__('snake.modules.fastqc.cmd_' + name,
                             None, None, ['cli'])
            
            return mod.cli
        return mod.cli

@click.group(cls=ComplexCLI, chain = True)
# @pass_environment
@click.pass_obj
def complex_gr(config):
    click.echo('complex_gr')
    pass



@cli.group(name='test')
def test_group():
    pass

test_group.add_command(complex_gr)

def main():
    print('======IN MAIN===========')
    print(argv[:])
    try:
        print(argv[:].index('run'))
    except:
        print('nnn')
    cli()