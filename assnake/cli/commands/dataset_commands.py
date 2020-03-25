import click, sys, os, glob, yaml, shutil
import pandas as pd
import assnake.api.loaders
import assnake.core.sample_set
from tabulate import tabulate
from assnake.api import fs_helpers
from assnake.utils import pathizer, dict_norm_print, download_from_url
from zipfile import ZipFile
from assnake.api.loaders import update_fs_samples_csv

# some util
def show_av_dict(dfs):
    '''
    print out available datasets from dict of db's (assnake.api.loaders.load_dfs_from_db(''))
    '''
    avail_dfs = ''
    for i, item in enumerate(dfs.keys()):
        avail_dfs += '{}. {} \t'.format(i + 1, item)

    if avail_dfs == '':
        avail_dfs = 'NA'
    click.echo('Available: ' + avail_dfs)
    exit(2)


# ---------------------------------------------------------------------------------------
#                                     LIST
# ---------------------------------------------------------------------------------------
@click.command(name='list')
def df_list():
    """List datasets in database"""
    dfs = assnake.api.loaders.load_dfs_from_db('')
    if len(list(dfs.keys())) == 0:
        click.echo('No datasets in your system yet!\nYou can create one by running\n' +
                   click.style('  assnake dataset create  ', bg='blue', fg='white', bold=True))

    for df in dfs.values():
        df_name = df['df']
        click.echo(click.style('' * 2 + df_name + ' ' * 2, fg='green', bold=True))
        click.echo('  Filesystem prefix: ' + df.get('fs_prefix', ''))
        click.echo('  Full path: ' + os.path.join(df.get('fs_prefix', ''), df['df']))
        click.echo('  Description: ')
        dict_norm_print(df.get('description', ''))
        click.echo('')


# ---------------------------------------------------------------------------------------
#                                   CREATE
# ---------------------------------------------------------------------------------------

def check_df_name_avialability(df):
    
    dfs_in_db = assnake.api.loaders.load_dfs_from_db('')
    if df in dfs_in_db:
        click.secho('Dataset with such name already exists!. Aborting.')
        exit()
    return df

def validate_df_prefix(df, fs_prefix, create_prefix = False, allow_not_empty_full_path = False):
    """
    Ensures that can write to full_df_path and it is empty
    """
    df = check_df_name_avialability(df)
    if df is None:
        exit()

    if not os.path.isabs(fs_prefix):
        click.secho('fs_prefix is not absolute!') 
        exit()
    full_df_path = os.path.join(fs_prefix, df)

    if os.path.isdir(full_df_path):
        if not os.listdir(full_df_path) and os.access(full_df_path, os.W_OK):
            click.echo('Great! ' + full_df_path + ' exists, is empty and assnake can write to that folder, initializing dataset here')
            return full_df_path
        elif os.listdir(full_df_path):
            if allow_not_empty_full_path:
                click.echo(full_df_path + ' exists, but it is not empty. allow_not_empty_full_path is true, existance logic will be handled downstream.')
                return full_df_path
            else:
                click.echo('Oops =( ' + full_df_path + ' exists, but it is not empty. Aborting.')
                return None
        elif not os.access(full_df_path, os.W_OK):
            click.echo('Oops =( ' + full_df_path + ' exists, but assnake doesnt have permissions to write there. Aborting.')
            return None
    else:
        click.echo('Oops =( ' + full_df_path + ' exists, but it is not directory. Aborting.')

        
    if os.path.isdir(fs_prefix):
        if os.access(fs_prefix, os.W_OK):
            click.echo('Great! ' + fs_prefix + ' exists, and assnake can write to that folder, initializing dataset in:')
            click.echo(full_df_path + '  (Will be created)')
            return full_df_path
    else:
        if not create_prefix:
            click.echo('Oops =( ' + fs_prefix + ' doesnt exist, and create_prefix is set to False. Aborting.')
        else:
            click.echo('Great! ' + fs_prefix + ' will be created, and dataset will be initialized in:')
            click.echo(full_df_path + '  (Will be created)')
            return full_df_path

    return None

def check_absolute_path(path, path_name):
    path = path.rstrip('\/')
    if not os.path.isabs(path):
        click.secho(path_name + ' is not absolute!') 
        exit()
    return path

@click.command(name='create')

@click.option('--df-name', '-d', help='Name of the dataset', required = False)

@click.option('--data-storage-folder', '-f', 
    help='Folder where you want to store your data. \
        Directory with provided df-name will be created inside this folder. MUST exist, may be not empty.\
        Resulting {storage-folder}/{df-name} directory MUST be empty, may not exist. \
        For regestiring existing datasets inside new Assnake instance use assnake dataset register', required=False)

@click.option('--full-path-to-df', '-p', 
    help='Full path to dataset folder. \
    Last part of the path (basename) will be used as df-name. MUST be empty, may not exist. \
    For regestiring existing datasets inside new Assnake instance use assnake dataset register', required=False)


@click.option('--first-preprocessing-name', 
    help='Name of your first preprocessing. raw by default. You want to set it to sra, for exaple if yoo are planning to download from NCBI. Purely cosmetic effect.', required=False, default = 'raw')

@click.option('--description', '-D', nargs=2, multiple=True, required=False, type=click.Tuple([str, str]),
              help='Add some description in this way ` assnake dataset create ... -D property_1 value_1 ... -D property_n value_n`')

@click.option('--quietly', '-q', is_flag=True, help='Doing it quietly. No questions.')
@click.option('--test-data', '-t', is_flag=True, help='Download test data from Humann2 tutorial')
@click.pass_obj
def df_create(config, df_name, data_storage_folder, full_path_to_df,first_preprocessing_name, description, quietly, test_data):
    """Register your dataset inside ASSNAKE!\n
        You can use it in interactive mode.
        Usage: assnake dataset create [dataset] or -d [dataset] ..

    """
    
    if ((df_name is None) or (data_storage_folder is None)):
        if full_path_to_df is None:
            click.echo('You must specify df_name AND data_storage_folder, OR full_path_to_df')
            exit()
        
        else: # full_path_to_df road
            # click.echo('df-name or data-storage-folder not set, using full_path_to_df')
            click.echo('Using full_path_to_df')
            full_path_to_df = check_absolute_path(full_path_to_df, 'full_path_to_df')
            df = os.path.basename(full_path_to_df)
            fs_prefix = os.path.dirname(full_path_to_df)
            full_df_path = validate_df_prefix(df, fs_prefix, True)

    # df_name AND data_storage_folder road
    elif (df_name is not None) and (data_storage_folder is not None):
        if full_path_to_df is not None:
            click.echo('df-name AND data-storage-folder are set, ignoring full_path_to_df')

        click.secho('Using df_name AND data_storage_folder')
        data_storage_folder = check_absolute_path(data_storage_folder, 'data_storage_folder')
        df = df_name # TODO VALIDATE
        fs_prefix = data_storage_folder
        full_df_path = validate_df_prefix(df, fs_prefix, True)

    os.makedirs(full_df_path, exist_ok=True)
    os.makedirs(os.path.join(full_df_path, 'reads', first_preprocessing_name), exist_ok=True)
    df_info = {'df': df, 'fs_prefix': fs_prefix, 'description': {}}

    assnake_db = config['config']['assnake_db']
    os.makedirs(os.path.join(assnake_db, 'datasets'), exist_ok=True)
    df_path_in_assnake = os.path.join(assnake_db, 'datasets', df)
    os.symlink(full_df_path, df_path_in_assnake, target_is_directory = True)

    with open(os.path.join(df_path_in_assnake, 'df_info.yaml'), 'w') as info_file:
        yaml.dump(df_info, info_file, default_flow_style=False)
    click.secho('Saved dataset ' + df + ' sucessfully!', fg='green')


    if test_data:
        download_from_url('http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/mgs_tutorial_Oct2017.zip', 
                                os.path.join(fs_prefix, df,'mgs_tutorial_Oct2017.zip'))
        
        with ZipFile(os.path.join(fs_prefix, df,'mgs_tutorial_Oct2017.zip'), 'r') as zipObj:
            
            zipObj.extractall(os.path.join(fs_prefix, df,'./'))
        shutil.rmtree(os.path.join(fs_prefix, df, 'reads'), ignore_errors=True)
        os.makedirs(os.path.join(fs_prefix, df, 'reads'), exist_ok=True)
        shutil.move (os.path.join(fs_prefix, df,'mgs_tutorial_Oct2017/raw_data'), os.path.join(fs_prefix, df, 'reads/raw'))
        # TODO cleanup
        shutil.rmtree(os.path.join(fs_prefix, df,'mgs_tutorial_Oct2017'), ignore_errors=True)

# ---------------------------------------------------------------------------------------
#                                   INIT
# ---------------------------------------------------------------------------------------

@click.command(name='init', help = 'Register dataset in Assnake based on the folder from where you called the command. (Working directory)')
@click.option('--df-name', '-d', help='Name of the dataset. If provided, folder with this name will be created in current dir.', required = False)
@click.option('--first-preprocessing-name', 
    help='Name of your first preprocessing. raw by default. You want to set it to sra, for exaple if yoo are planning to download from NCBI. Purely cosmetic effect.', required=False, default = 'raw')

@click.pass_obj
def df_init(config, df_name, first_preprocessing_name):
    cwd = os.getcwd()

    if df_name is None:
        df = os.path.basename(cwd)
        fs_prefix = os.path.dirname(cwd)
        full_df_path = validate_df_prefix(df, fs_prefix, True, True)
    else:
        df = df_name
        fs_prefix = cwd
        full_df_path = validate_df_prefix(df, fs_prefix, True, True)

    full_path_empty = True
    if os.path.isdir(full_df_path) and not os.listdir(full_df_path):
        click.echo(full_df_path + ' not empty!')
        click.echo('Trying to import as existing Assnake dataset (Not properly implemented)')
        full_path_empty = False
 

    os.makedirs(os.path.join(full_df_path, 'reads', first_preprocessing_name), exist_ok=True)
    df_info = {'df': df, 'fs_prefix': fs_prefix, 'description': {}}

    assnake_db = config['config']['assnake_db']
    os.makedirs(os.path.join(assnake_db, 'datasets'), exist_ok=True)
    df_path_in_assnake = os.path.join(assnake_db, 'datasets', df)
    os.symlink(full_df_path, df_path_in_assnake, target_is_directory = True)

    df_info_loc = os.path.join(df_path_in_assnake, 'df_info.yaml')
    if not os.path.isfile(df_info_loc):
        with open(df_info_loc, 'w') as info_file:
            yaml.dump(df_info, info_file, default_flow_style=False)
    click.secho('Saved dataset ' + df + ' sucessfully!', fg='green')

# ---------------------------------------------------------------------------------------
#                                   INFO
# ---------------------------------------------------------------------------------------
@click.command(name='info')
@click.option('--df', '-d', help='Name of the dataset', required=False)
@click.option('--preproc', '-p', help='Show samples for preprocessing', required=False)
@click.argument('df_arg', required=True)
@click.pass_obj
def df_info(config, df, preproc, df_arg):
    """View info for the specific dataset
        Usage: assnake dataset info [dataset] or -d [dataset] ...

    """
    if not (bool(df is None) ^ bool(df_arg is None)):
        click.echo('Please, specify dataset either as option or argument')
        df = click.prompt('Type the name in')
    if df is None:
        df = df_arg
    dfs = assnake.api.loaders.load_dfs_from_db('')
    try:
        df_info = dfs[df]
    except KeyError as e:
        click.echo('Can`t reach database with such name')
        show_av_dict(dfs)
    click.echo(click.style('' * 2 + df_info['df'] + ' ' * 2, fg='green', bold=True))
    click.echo('Filesystem prefix: ' + df_info.get('fs_prefix', ''))
    click.echo('Full path: ' + os.path.join(df_info.get('fs_prefix', ''), df))
    click.echo('Description: ')
    dict_norm_print(df_info.get('description', ''))

    mg_samples_loc = os.path.join(config['config']['assnake_db'], 'datasets', df_info['df'], 'mg_samples.tsv')

    if os.path.isfile(mg_samples_loc):
        click.echo(tabulate(pd.read_csv(mg_samples_loc, sep='\t'), headers='keys', tablefmt='fancy_grid'))

    reads_dir = os.path.join(df_info['fs_prefix'], df_info['df'], 'reads/*')
    preprocs = [p.split('/')[-1] for p in glob.glob(reads_dir)]
    preprocs.sort()
    preprocessing = {}
    all_samples = []
    for p in preprocs:
        samples = assnake.core.sample_set.SampleSet(df_info['fs_prefix'], df_info['df'], p)
        # samples.add_samples(df_info['fs_prefix'], df_info['df'], p)
        all_samples += (list(samples.samples_pd['fs_name']))
        samples = samples.samples_pd[['fs_name', 'reads']].to_dict(orient='records')
        # click.secho(p + ' ' + str(len(samples)) + ' samples')
        preprocessing.update({p: samples})

    click.echo('\nTotal samples: ' + str(len(set(all_samples))) + '\n')

    for key, value in preprocessing.items():
        click.echo('Samples in ' + click.style(key, bold=True) + ': ' + str(len(value)))
    if preproc is not None:
        samples_pd = pd.DataFrame(preprocessing[preproc])
        # print(samples_pd)
        click.echo(tabulate(samples_pd.sort_values('reads'), headers='keys', tablefmt='fancy_grid'))


# ---------------------------------------------------------------------------------------
#                                   DELETE
# ---------------------------------------------------------------------------------------
@click.command(name='delete')
@click.option('--df', '-d',
              help='Name of dataset to delete.', type=click.STRING )
@click.option('--hard', help='If is set, hard removing will be used instead of modifying congig file', is_flag=True)
@click.argument('df_arg', required=False)
@click.pass_obj
def df_delete(config, df, hard, df_arg):
    """
    Delete datasets
    Usage: assnake dataset delete [dataset] or -d [dataset] ...
    """
    if not (bool(df is None) ^ bool(df_arg is None)):
        click.echo('Please, specify dataset either as option or argument')
        df = click.prompt('Type the name in:')
    if df is None:
        df = df_arg
    dfs = assnake.api.loaders.load_dfs_from_db('')
    try:
        df_info = dfs[df]
    except KeyError as e:
        click.echo('Can`t reach database with such name')
        show_av_dict(dfs)

    respond = fs_helpers.delete_ds(df)
    if respond[0]:
        click.secho('Successfully deleted', fg='bright_white', bg='green')
    else:
        click.secho('ERROR', bg='red')
        click.echo('For details see traceback below')
        click.echo(respond[1])

    if hard and click.confirm(
        'Are you sure to delete this nice and probably huge datasets, which you may redownload for eternity -- use modifying config instead?',
        abort=True):
        shutil.rmtree(os.path.join(df_info.get('fs_prefix', ''), df))

# ---------------------------------------------------------------------------------------
#                                   IMPORT-READS
# ---------------------------------------------------------------------------------------
# DONE decide if we need either d and t or proceed both arguments as one and automatically choose path or not
@click.command(name='import-reads')
@click.option('--reads', '-r', prompt='Location of folder with read files',
              help='Location of folder with read files', type=click.Path())
@click.option('--dataset', '-d', help='Assnake dataset name. If -t is not specified', required=False)
@click.option('--rename-method', help='How to rename samples', type=click.Choice(['replace-', 'removeSending'], case_sensitive=False), required=False)
@click.option('--target', '-t', help='Location of the target directory. If -d is not specified.', required=False,
              type=click.Path())
@click.option('--sample_set', '-s', help='Comma-divided list of samples of interest', required=False)
@click.option('--sample_list', '-l', help='Location of file with line by line samples of interest', required=False,
              type=click.Path())
@click.option('--copy', help='If is set, hard copying will be used instead of symbolic links ', is_flag=True)
@click.pass_obj
def df_import_reads(config, reads, dataset, rename_method, target, sample_set, sample_list, copy):
    """
    Import reads from directory to assnake dataset. Currently local text files are supported. The --target argument
    point to location (relative or absolute) of assnake dataset in your file system. Please, pay attention,
    that -t and -d  arguments are excclusive for each over -- specify only one of them -- as well as -s and -l.
    With -s `sample_1,sample_2,...,sample_n` notation is implied (no whitespaces between sample names)
    """
    # print(os.getcwd())
    # This is cli wrapping over create_links from fs_helpers

    # stuff about presence of arguments
    arg_d = not bool(dataset is None)
    arg_t = not bool(target is None)
    arg_s = not bool(sample_set is None)
    arg_l = not bool(sample_list is None)

    # check if samples arguments are ok
    if arg_l & arg_s:
        click.secho('Collision tends to be observed. Please, specify either list of samples in prompt or in file',
                    err=True)
        exit(1)

    # check if destination args are ok
    if not (arg_d ^ arg_t):
        click.secho('Please, specify either database (-d) or absolute path (-t)', err=True)
        exit(1)

    # some stuffff to ensure correctness of source and destination (how philosophical)
    if arg_d:
        try:
            df_info = assnake.api.loaders.load_df_from_db(dataset)
        except Exception as e:
            dfs = assnake.api.loaders.load_dfs_from_db('')
            click.echo('Can`t reach database with such name', err=True)
            show_av_dict(dfs)
        target = '{}/{}/reads/raw'.format(df_info['fs_prefix'], df_info['df'])
    else:
        target = pathizer(target)
        if not os.path.exists(target):
            click.secho("Provided sample-list file couldn't be detected", err=True)
            exit(2)
    target = '{}/{}/reads/raw'.format(df_info['fs_prefix'], df_info['df'])
    os.makedirs(target, exist_ok=True)
    # def rename(sample):
    #     new_name_wc = 'Rst{}_{}_B2'
    #     splitted = sample.split('_')
    #     splitted[2] = (1 if len(splitted) == 3 else a)
    #     return new_name_wc.format(splitted[1], splitted[2])

    if rename_method == 'removeSending':
        modify_name=lambda arg: '_'.join(arg.replace('-', '_').split('_')[0:-1])
    else:
        modify_name=lambda arg: arg.replace('-', '_')

    samples_in_run = fs_helpers.get_samples_from_dir(reads, modify_name)
    samples_in_run['fs_name'] = samples_in_run['modified_name']
    fs_helpers.create_links(target,  samples_in_run, hard=copy)

    update_fs_samples_csv(df_info['df'])
    click.secho("SUCCESSFULLY IMPORTED READS!", bg='green') 
    # sample_names = {d['sample_name'] for d in dicts}
    # if arg_s:
    #     samples_of_interest = sample_set.split(',')
    # elif arg_l:
    #     sample_list = pathizer(sample_list)
    #     if os.path.exists(sample_list):
    #         with open(sample_list, 'r') as ls_file:
    #             samples_of_interest = ls_file.readlines()
    #     else:
    #         click.secho("Provided sample-list file couldn't be detected", err=True)
    # else:
    #     samples_of_interest = sample_names

    # # Such a mess! check if specified samples are in directory
    # if (arg_s or arg_l) and (len(set(samples_of_interest) - {d['sample_name'] for d in dicts}) != 0):
    #     click.secho("Warning! Samples, been specified, are not in reads file", err=True)
    #     click.confirm('Continue?', abort=True)
    # if copy:
    #     click.secho("Please, keep in mind, that hard copying may take plenty of time")

    # Here the logic -- just call create_links
    


@click.command(name='rescan')
@click.option('--dataset', '-d', help='Assnake dataset name', required=False)
@click.argument('df_arg', required=False)
@click.pass_obj
def rescan_dataset(config, dataset, df_arg):
    """
    Now it just updates fs_samples.tsv in ./assnkae_db/{dataset}/

    Usage: assnake dataset rescan [dataset] or -d [dataset] ..
    """
    if not (bool(dataset is None) ^ bool(df_arg is None)):
        click.echo('Please, specify dataset either as option or argument')
        dataset = click.prompt('Type the name in:')
    if dataset is None:
        dataset = df_arg
    success = update_fs_samples_csv(dataset)
    if success:
        click.secho('SUCCESSFULLY UPDATED INFORMATION IN DATABASE!', fg='green')
