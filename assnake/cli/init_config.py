import os, yaml
from assnake.utils import read_yaml


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
