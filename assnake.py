import click, sys, os, glob, yaml
import pandas as pd
import api.loaders
import api.sample_set
from tabulate import tabulate
import snakemake

@click.group()
@click.version_option()
@click.pass_context
def cli(ctx):
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_loc = os.path.join(dir_of_this_file, 'config.yml')
    config = {}
    with open(config_loc, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    wc_config_loc = os.path.join(dir_of_this_file, 'wc_config.yaml')
    wc_config = {}
    with open(wc_config_loc, 'r') as stream:
        try:
            wc_config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    ctx.obj = {'config': config, 'wc_config': wc_config}
    pass #Entry Point

@cli.group()
def dataset():
    """Commands to work with datasets"""
    pass

@click.command()
def df_list():
    """List datasets in database"""
    dfs = api.loaders.load_dfs_from_db('')
    for df in dfs.values():
        df_name = df['df']
        click.echo(click.style(''*2 + df_name + ' '*2, fg='green', bold=True))
        click.echo('  Filesystem prefix: ' + df.get('fs_prefix', ''))
        click.echo('  Description: ' + df.get('description', ''))
        click.echo('')

@click.command()
@click.option('--df','-d', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--fs_prefix','-p', prompt='Filesystem prefix', help='Filesystem prefix' )
@click.pass_obj
def create_df(config, df, fs_prefix):
    """Create entry for dataset in database."""
    assnake_db_search = os.path.join(config['config']['assnake_db'], 'datasets/*')
    dfs = [d.split('/')[-1] for d in glob.glob(os.path.join(assnake_db_search))]
    print(df)
    print(fs_prefix)
    if df not in dfs:
        if os.path.isdir(os.path.join(fs_prefix, df)):
            df_info = {'df': df, 'fs_prefix': fs_prefix}
            os.mkdir(os.path.join(config['config']['assnake_db'], 'datasets/'+df))
            with open(os.path.join(config['config']['assnake_db'], 'datasets/'+df,'df_info.yaml'), 'w') as info_file:
                yaml.dump(df_info, info_file, default_flow_style=False)
            click.secho('Saved dataset ' + df + ' sucessfully!', fg='green')
    else:
        click.secho('Duplicate name!', fg='red')

@click.command()
@click.option('--name','-n', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--preproc','-p', help='Show samples for preprocessing', required=False)
def df_info(name, preproc):
    """View info for the specific dataset"""
    print(preproc)

    dfs = api.loaders.load_dfs_from_db('')
    df_info = dfs[name]
    click.echo(click.style(''*2 + df_info['df'] + ' '*2, fg='green', bold=True))
    click.echo('Filesystem prefix: ' + df_info.get('fs_prefix', ''))
    click.echo('Full path: ' + os.path.join(df_info.get('fs_prefix', ''), name))
    click.echo('Description: ' + df_info.get('description', ''))

    reads_dir = os.path.join(df_info['fs_prefix'], df_info['df'], 'reads/*')
    preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
    preprocs.sort()
    preprocessing = {}
    all_samples = []
    for p in preprocs:
        samples = api.sample_set.SampleSet()
        samples.add_samples(df_info['fs_prefix'], df_info['df'], p)
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
dataset.add_command(create_df)

@cli.group()
def result():
    """Commands to work with results"""
    pass

@click.command()
@click.option('--df','-d', prompt='Name of the dataset', help='Name of the dataset' )
@click.option('--preproc','-p', prompt='Preprocessing to use', help='Preprocessing to use' )
@click.option('--result','-r', prompt='Result you want', help='Preprocessing to use' )

@click.option('--params', help='Parameters id to use', required=False, default = 'def')
@click.option('--list-name', help='Name of sample set if required by rule')

@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--run/--no-run', default=False)
@click.pass_obj
def request(config, df, preproc, result, params, list_name,   threads, jobs, run):
    """Request result for your samples.\nYou need to provide info defining sample set, or run it from {prefix}/{df}/reads/{preproc} directory"""

    ss = api.sample_set.SampleSet()
    df = api.loaders.load_df_from_db(df)
    print(df)
    ss.add_samples(df['fs_prefix'], df['df'], preproc)
    click.echo(tabulate(ss.samples_pd[['fs_name', 'reads', 'preproc']].sort_values('reads'), headers='keys', tablefmt='fancy_grid'))

    res_list = []

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

    elif result=='mp2':
        res_list = ss.get_locs_for_result(result)
    elif result=='fastqc':
        res_list = ss.get_locs_for_result(result)
    elif result == 'dada2-filter-and-trim':
        res_list = ss.get_locs_for_result(result, params=params)
    elif result=='dada2-learn-errors':
        # ss.prepare_dada2_sample_list()
        if list_name is None:
            list_name = click.prompt('Please enter name for your dada2 error profile')
            print(list_name)
        ss.prepare_dada2_sample_list(list_name)
        res_list= [(config['wc_config']['dada2_err_wc'].format(
                dada2_dir = config['config']['dada2_dir'],
                sample_set = list_name,
                strand = strand,
                params=params
            )) for strand in ['R1', 'R2']]
        

    
    curr_dir = os.path.abspath(os.path.dirname(__file__))
    click.echo(curr_dir)
    click.echo(jobs)

    status = snakemake.snakemake(os.path.join(curr_dir, './bin/snake/base.py'), 
        config = dict(assnake_install_dir= os.path.abspath(os.path.dirname(__file__))),    
        targets=res_list, 
        printshellcmds=True,
        dryrun=not run, 
        configfile=os.path.join(curr_dir, 'config.yml'),
        drmaa_log_dir = os.path.join(curr_dir, '../pipeline_dev/drmaa_log'),
        use_conda = True,
        conda_prefix = os.path.join(curr_dir, '../pipeline_dev/conda'),
        drmaa=' -V -S /bin/bash -pe make {threads}'.format(threads=threads),
        cores=jobs, nodes=jobs
        )

        

result.add_command(request)

if __name__ == '__main__':
    cli()