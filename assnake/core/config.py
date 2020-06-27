import yaml, configparser, os, click
import os, sys
import requests, urllib
from tqdm import tqdm
from assnake.utils.general import read_yaml
from pathlib import Path

def load_wc_config():
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    return (read_yaml(os.path.join(dir_of_this_file, '../snake/wc_config.yaml')))


def read_internal_config():
    '''
    Reads the config at ~/.config/assnake/internal_config.yaml and returns as dict
    '''
    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')
    
    if os.path.isfile(internal_config_loc):
        internal_config = read_yaml(internal_config_loc)
        return internal_config
    else:
        internal_config_dir = os.path.join(str(Path.home()), '.config/assnake/')

        os.makedirs(internal_config_dir, exist_ok=True)

        if not os.path.isfile(internal_config_loc):
            with open(internal_config_loc, 'w+') as file:
                _ = yaml.dump({'instance_config_loc': 'not_set'}, file, sort_keys=False)
        return {'instance_config_loc': ''}

def read_assnake_instance_config():
    '''
    Reads particular assnake instance config. It is stored inside assnake database as config.yaml (Name subject to change). 
    :return: Returns dict if instance config exists, None otherwise.
    '''
    internal_config = read_internal_config()
    instance_config_loc = internal_config['instance_config_loc']

    if os.path.isfile(instance_config_loc):
        return read_yaml(instance_config_loc)
    else:
        return None

def update_internal_config(update_dict):
    '''
    Updates internal config with provided dictionary

    :return: updated internal config if everything is ok, None otherwise
    '''
    internal_config = read_internal_config()
    internal_config.update(update_dict)
    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')
    with open(internal_config_loc, 'w+') as file:
        _ = yaml.dump(internal_config, file, sort_keys=False)

    return internal_config

def check_if_assnake_is_initialized():
    '''
    Just tries to read assnake instance config, and prints an error if there is no instance config.
    '''
    instance_config = read_assnake_instance_config()

    if instance_config is None:
        click.secho("You need to init your installation!", fg='red', bold=True)
        click.echo("Don't worry, it won't take long.")
        click.echo('Just run ' + click.style('assnake init start', bg='blue', fg='bright_white'))
        exit()




## ====RENAME THIS === update_instance_config
def update_config(update_dict):
    '''
    Updates instance config with provided dictionary

    :return: updated instance config if everything is ok, None otherwise
    '''

    instance_config = read_assnake_instance_config()
    instance_config.update(update_dict)

    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')

    with open(internal_config_loc, 'w+') as file:
        _ = yaml.dump(instance_config, file, sort_keys=False)

    return instance_config


## ============================================

def fill_and_write_instance_config(assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, instance_config_location):
    '''
    Accepts all the necessary config props and dumps them to the provided location as file

    :return: location of instance config file if everything is ok, None otherwise
    '''
    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    
    assnake_install_dir = (dir_of_this_file.replace('assnake/commands', 'snake'))

    config_template_loc = os.path.join(dir_of_this_file, '../snake/config_template.yml')
    config_template = read_yaml(config_template_loc)

    config_template['assnake_db']          = assnake_db
    config_template['assnake_install_dir'] = assnake_install_dir
    config_template['fna_db_dir']          = fna_db_dir
    config_template['bwa_index_dir']       = bwa_index_dir
    config_template['conda_dir']           = conda_dir
    config_template['drmaa_log_dir']       = drmaa_log_dir

    os.makedirs(fna_db_dir, exist_ok=True)
    os.makedirs(bwa_index_dir, exist_ok=True)
    os.makedirs(conda_dir, exist_ok=True)
    os.makedirs(drmaa_log_dir, exist_ok=True)

    with open(instance_config_location, 'w+') as file:
        _ = yaml.dump(config_template, file, sort_keys=False)

    return instance_config_location