import click, os, snakemake
# from assnake.assnake_cli import pass_environment


@click.command('snake', short_help='Runs snakemake for requested results')

@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--run/--no-run', default=False)

@click.pass_obj

def cli(config, threads, jobs, run):
    
    print(config['requests'])
    click.echo('RUN SNAKEMAKE')

    curr_dir = os.path.abspath(os.path.dirname(__file__))
    # click.echo('Current dir: ' + curr_dir)
    # click.echo('Number of jobs torun in parallel: ' + str(jobs))

    status = snakemake.snakemake(os.path.join(curr_dir, '../bin/snake/base.py'), 
        # config = config['wc_config'],    
        targets=config['requests'], 
        printshellcmds=True,
        dryrun=not run, 
        configfiles=[os.path.join(curr_dir, '../../snake/config.yml')],
        drmaa_log_dir = config['config']['drmaa_log_dir'],
        use_conda = True,
        latency_wait = 120,
        conda_prefix = config['config']['conda_dir'],
        drmaa=' -V -S /bin/bash -pe make {threads}'.format(threads=threads),
        cores=jobs, nodes=jobs
        )