import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.api.sample_set
from tabulate import tabulate
import snakemake

@click.group()
@click.version_option()
@click.pass_context
def cli(ctx):
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_loc = os.path.join(dir_of_this_file, '../snakemake/config.yml')

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

    wc_config_loc = os.path.join(dir_of_this_file, '../snakemake/wc_config.yaml')
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
    config_loc = os.path.join(dir_of_this_file, '../snakemake/config.yml')

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

@cli.group()
def dataset():
    """Commands to work with datasets"""
    pass

@click.command(name='list')
def df_list():
    """List datasets in database"""
    dfs = assnake.api.loaders.load_dfs_from_db('')
    if len(list(dfs.keys())) == 0:
        click.echo('No datasets in your system yet!\nYou can create one by running\n' + 
            click.style('  assnake dataset create  ', bg='blue', fg='white', bold=True))

    for df in dfs.values():
        df_name = df['df']
        click.echo(click.style(''*2 + df_name + ' '*2, fg='green', bold=True))
        click.echo('  Filesystem prefix: ' + df.get('fs_prefix', ''))
        click.echo('  Full path: ' + os.path.join(df.get('fs_prefix', ''), df['df']))
        click.echo('  Description: ' + df.get('description', ''))
        click.echo('')

@click.command(name='create')
@click.option('--df','-d', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--fs_prefix','-p', prompt='Filesystem prefix', help='Filesystem prefix' )
@click.pass_obj
def df_create(config, df, fs_prefix):
    """Create entry for dataset in database."""
    assnake_db_search = os.path.join(config['config']['assnake_db'], 'datasets/*')
    dfs = [d.split('/')[-1] for d in glob.glob(os.path.join(assnake_db_search))]
    print(df)
    print(fs_prefix)
    if df not in dfs:
        if os.path.isdir(os.path.join(fs_prefix, df)):
            df_info = {'df': df, 'fs_prefix': fs_prefix}
            os.makedirs(os.path.join(config['config']['assnake_db'], 'datasets/'+df), exist_ok=True)
            with open(os.path.join(config['config']['assnake_db'], 'datasets/'+df,'df_info.yaml'), 'w') as info_file:
                yaml.dump(df_info, info_file, default_flow_style=False)
            click.secho('Saved dataset ' + df + ' sucessfully!', fg='green')
    else:
        click.secho('Duplicate name!', fg='red')

@click.command(name ='info')
@click.option('--name','-n', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--preproc','-p', help='Show samples for preprocessing', required=False)
@click.pass_obj
def df_info(config, name, preproc):
    """View info for the specific dataset"""

    dfs = assnake.api.loaders.load_dfs_from_db('')
    df_info = dfs[name]
    click.echo(click.style(''*2 + df_info['df'] + ' '*2, fg='green', bold=True))
    click.echo('Filesystem prefix: ' + df_info.get('fs_prefix', ''))
    click.echo('Full path: ' + os.path.join(df_info.get('fs_prefix', ''), name))
    click.echo('Description: ' + df_info.get('description', ''))
    
    mg_samples_loc = os.path.join(config['config']['assnake_db'], 'datasets', df_info['df'], 'mg_samples.tsv')

    if os.path.isfile(mg_samples_loc):
        click.echo(tabulate(pd.read_csv(mg_samples_loc, sep='\t'), headers='keys', tablefmt='fancy_grid'))

    reads_dir = os.path.join(df_info['fs_prefix'], df_info['df'], 'reads/*')
    preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
    preprocs.sort()
    preprocessing = {}
    all_samples = []
    for p in preprocs:
        samples = assnake.api.sample_set.SampleSet(df_info['fs_prefix'], df_info['df'], p)
        # samples.add_samples(df_info['fs_prefix'], df_info['df'], p)
        all_samples += (list(samples.samples_pd['fs_name']))
        samples = samples.samples_pd[['fs_name', 'reads']].to_dict(orient='records')
        # click.secho(p + ' ' + str(len(samples)) + ' samples')
        preprocessing.update({p:samples})

    click.echo('\nTotal samples: ' + str(len(set(all_samples))) + '\n')


    for key, value in preprocessing.items():
        click.echo('Samples in ' + click.style(key, bold=True) + ': ' + str(len(value)))
    if preproc is not None:
        samples_pd = pd.DataFrame(preprocessing[preproc])
        # print(samples_pd)
        click.echo(tabulate(samples_pd.sort_values('reads'), headers='keys', tablefmt='fancy_grid'))


dataset.add_command(df_list)
dataset.add_command(df_info)
dataset.add_command(df_create)

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
    print(samples_to_add)
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
            while click.confirm('Do you want to process any sets?'):
                value = click.prompt('Please, enter set number', type=int)
                if value+1 > len(prepared_sets):
                    click.echo(click.style('NO SUCH SAMPLE SET',fg='red'))
                else:
                    selected_sets.append(value)
                    click.echo(click.style('Selected sets: ' + str(set(selected_sets)), fg='green'))
            min_len = click.prompt('Please, enter minimum contig lenth', type=int)

            for ss in selected_sets:
                res_list += [prepared_sets[ss].replace('sample_set.tsv', 'final_contigs__{min_len}.fa'.format(min_len=min_len))]
        elif result=='metabat2':
            mb2 = '{fs_prefix}/{df}/metabat2/bwa__0.7.17__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/metabat2.done'

            # prepared_sets = glob.glob(os.path.join(df['fs_prefix'],df['df'],'assembly/*/*/sample_set.tsv'))
            # click.echo('We found ' + str(len(prepared_sets)) + ' prepared sets:')
            # for i, pset in enumerate(prepared_sets):
            #     name = pset.split('/')[-2]
            #     run_info = pset.split('/')[-3]
            #     click.echo(click.style(str(i), bold = True) + ' ' + run_info + ' ' + name)

            # selected_sets = []
            # while click.confirm('Do you want to process any sets?'):
            #     value = click.prompt('Please, enter set number', type=int)
            #     if value+1 > len(prepared_sets):
            #         click.echo(click.style('NO SUCH SAMPLE SET',fg='red'))
            #     else:
            #         selected_sets.append(value)
            #         click.echo(click.style('Selected sets: ' + str(set(selected_sets)), fg='green'))
            # min_len = click.prompt('Please, enter minimum contig lenth', type=int)

            # for ss in selected_sets:
            #     name = pset.split('/')[-2]
            #     name = 'p136'
            #     res_list += [mb2.format(fs_prefix=df['fs_prefix'],df=df['df'], params='def',sample_set=name, mod=1000)]
            name = 'p136'
            res_list += [mb2.format(fs_prefix=df['fs_prefix'],df=df['df'], params='def',sample_set=name, mod=1000)]
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
        configfiles=[os.path.join(curr_dir, '../snakemake/config.yml')],
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
                print(filename)
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        try:
            if sys.version_info[0] == 2:
                name = name.encode('ascii', 'replace')
            print('assnake.commands.cmd_' + name)
            mod = __import__('assnake.commands.cmd_' + name,
                             None, None, ['cli'])
        except ImportError as error:
            print('import_error')
            print(error.__class__.__name__ + ": ")
            return
        return mod.cli

@click.command(cls=ComplexCLI)
# @pass_environment
@click.pass_obj
def complex_gr(config):
    click.echo('results')
    pass



@cli.group(name='test')
def test_group():
    pass

result.add_command(complex_gr)


if __name__ == '__main__':
    cli()