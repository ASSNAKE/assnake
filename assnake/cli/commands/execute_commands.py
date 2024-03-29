#############################################
#   Perform assnake  ***  RUN cli command   #
#############################################

import click, os
from assnake.core.config import read_internal_config


#---------------------------------------------------------------------------------------
#                                     RUN
#---------------------------------------------------------------------------------------

@click.command('gather', short_help='Gathers requested results and passes them to Snakemake')
@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--drmaa/--no-drmaa', default=False)
@click.option('--run/--no-run', default=False)
@click.option('--touch/--no-touch', default=False)
@click.option('--unlock/--no-unlock', default=False)
@click.option('--debug-dag/--no-debug-dag', default=False)
@click.option('--dag/--no-dag', default=False)
@click.option('--lint/--no-lint', default=False)
@click.option('--print-compilation/--no-print-compilation', default=False)
@click.pass_obj

def gather(config, threads, jobs, drmaa, run, touch, unlock, debug_dag, dag, lint, print_compilation):
    import snakemake # Moved import here because it is slow as fucking fuck
    
    
    internal_config = read_internal_config()
    # print(config['requests'])

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
        configfiles=[internal_config['instance_config_loc']],
        drmaa_log_dir = config['config']['drmaa_log_dir'],
        
        lint = lint,
        listrules = print_compilation,
        printdag   = dag, 

        conda_frontend='conda',
        use_conda = True,
        rerun_triggers = ['mtime'],
        conda_prefix = config['config']['conda_dir'],
        keepgoing = True,
        latency_wait = 10,
        drmaa=drmaa_param,
        touch = touch,
        unlock = unlock,
        debug_dag = debug_dag, 
        cores=jobs, nodes=jobs)

    # print(config['requested_results']) 
    
    if run:
        click.echo('Updating Datasets:' + str(config['requested_dfs']))
        # for requested_df in set(config['requested_dfs']):
        #     update_fs_samples_csv(requested_df)

        print(config['requested_results']) 
