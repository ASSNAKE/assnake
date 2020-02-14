#############################################
#   Perform assnake  ***  RUN cli command   #
#############################################

import click, os, snakemake
from assnake.utils import get_config_loc, load_config_file
from assnake.api.update_fs_samples import update_fs_samples_csv
# from assnake.assnake_cli import pass_environment


#---------------------------------------------------------------------------------------
#                                     RUN
#---------------------------------------------------------------------------------------

@click.command('run', short_help='Runs snakemake for requested results')
@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--drmaa/--no-drmaa', default=False)
@click.option('--run/--no-run', default=False)

@click.pass_obj

def run(config, threads, jobs, drmaa, run):
    print(config['requests'])

    click.secho('-----===RUN SNAKEMAKE===-----', bg='green', fg='black')

    curr_dir = os.path.abspath(os.path.dirname(__file__))
    # click.echo('Current dir: ' + curr_dir)
    # click.echo('Number of jobs torun in parallel: ' + str(jobs))
    drmaa_param = None
    if drmaa:
        drmaa_param=' -V -S /bin/bash -pe make {threads}'.format(threads=threads)
    status = snakemake.snakemake(os.path.join(curr_dir, '../../snake/snake_base.py'), 
        # config = config['wc_config'],    
        targets=config['requests'], 
        printshellcmds=True,
        dryrun=not run, 
        # config = load_config_file(),
        configfiles=[get_config_loc()],
        drmaa_log_dir = config['config']['drmaa_log_dir'],
        use_conda = True,
        latency_wait = 120,
        conda_prefix = config['config']['conda_dir'],
        drmaa=drmaa_param,
        cores=jobs, nodes=jobs)
    if run:
        update_fs_samples_csv(**config['invoked_args'])


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\