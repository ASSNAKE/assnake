#############################################
#   Perform assnake  ***  RUN cli command   #
#############################################

import click, os, yaml
from assnake.core.config import read_internal_config
from assnake.api.loaders import update_fs_samples_csv
from assnake.cli.commands.dataset_commands import df_update_info

#---------------------------------------------------------------------------------------
#                                     RUN
#---------------------------------------------------------------------------------------

@click.command('gather', short_help='Gathers requested results and passes them to Snakemake')
@click.option('--threads','-t', help='Threads per job', default=4)
@click.option('--jobs','-j', help='Number of jobs', default=1)
@click.option('--drmaa/--no-drmaa', default=False)
@click.option('--run/--no-run', default=False)
@click.option('--touch/--no-touch', default=False)
@click.pass_obj

def gather(config, threads, jobs, drmaa, run, touch):
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
        use_conda = True,
        latency_wait = 120,
        conda_prefix = config['config']['conda_dir'],
        drmaa=drmaa_param,
        touch = touch,
        cores=jobs, nodes=jobs,
        # ignore_ambiguity=True
        )
    # print('STATUS')
    # print(status)
    # print()
    # print(config['requests'])
    # print(config['requests_storage'])
    # print(config['requested_results']) 
    
    if run:
        click.echo('Updating Datasets:' + str(config['requested_dfs']))
        for requested_df in set(config['requested_dfs']):
            update_fs_samples_csv(requested_df)

        if 'requests_storage' in config:
            for df in config['requests_storage'].keys():
                df_info_loc = config['config']['assnake_db']+'/datasets/{df}/df_info.yaml'.format(df = df)
                df_info = {}
                if os.path.isfile(df_info_loc):
                    with open(df_info_loc, 'r') as stream:
                        try:
                            info = yaml.load(stream, Loader=yaml.FullLoader)
                            if 'df' in info:
                                df_info =  info
                        except yaml.YAMLError as exc:
                            pass
                if df_info is not {}:
                    processes = df_info['processes'] if 'processes' in df_info else {}
                    for process_name in config['requests_storage'][df]:
                        if process_name not in processes:
                            processes[process_name] = []
                        file_set = [filepath.replace(df_info['full_path']+'/', '') for filepath in config['requests_storage'][df][process_name]
                                    if os.path.exists(filepath)]
                        processes[process_name] += file_set
                    df_update_info(df_name = df, processes=processes, clean_paths=True)
        else:
            click.secho("df_info.yaml could not be update, since this feature was not implemented for this function.",  fg='yellow')
    else:
        if 'destroy_if_not_run' in config:
            if 'files' in config['destroy_if_not_run']:
                for file in config['destroy_if_not_run']['files']:
                    os.remove(file)
            if 'directories' in config['destroy_if_not_run']:
                for folder in config['destroy_if_not_run']['directories']:
                    if len(os.listdir(folder)) == 0:
                        os.rmdir(folder)


        # print(config['requests'])
        # print(config['requests_storage'])
        # print(config['requested_results']) 