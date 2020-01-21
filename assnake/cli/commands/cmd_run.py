import click, os, snakemake
from assnake.utils import get_config_loc
# from assnake.assnake_cli import pass_environment


@click.command('run', short_help='Runs snakemake for requested results')

@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--drmaa/--no-drmaa', default=False)
@click.option('--run/--no-run', default=False)

@click.pass_obj

def run(config, threads, jobs, drmaa, run):
    
    click.echo('RUN SNAKEMAKE')

    curr_dir = os.path.abspath(os.path.dirname(__file__))
    # click.echo('Current dir: ' + curr_dir)
    # click.echo('Number of jobs torun in parallel: ' + str(jobs))
    drmaa_param = None
    if drmaa:
        drmaa_param=' -V -S /bin/bash -pe make {threads}'.format(threads=threads)

    status = snakemake.snakemake(os.path.join(curr_dir, '../bin/snake/base.py'), 
        # config = config['wc_config'],    
        targets=config['requests'], 
        printshellcmds=True,
        dryrun=not run, 
        configfiles=[get_config_loc()],
        drmaa_log_dir = config['config']['drmaa_log_dir'],
        use_conda = True,
        latency_wait = 120,
        conda_prefix = config['config']['conda_dir'],
        drmaa=drmaa_param,
        cores=jobs, nodes=jobs)