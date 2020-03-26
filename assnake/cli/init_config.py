import os, yaml
from assnake.utils import read_yaml, get_internal_config


# @graph_of_calls('fill_n_write_config.png')
def fill_and_write_config(assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, config_location):
    # config_internal = get_internal_config()

    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    assnake_install_dir = (dir_of_this_file.replace('assnake/commands', 'snake'))

    config_template_loc = os.path.join(dir_of_this_file, '../snake/config_template.yml')

    config_template = read_yaml(config_template_loc)

    config_template['assnake_db'] = assnake_db
    config_template['assnake_install_dir'] = assnake_install_dir
    config_template['fna_db_dir'] = fna_db_dir
    config_template['bwa_index_dir'] = bwa_index_dir
    config_template['conda_dir'] = conda_dir
    config_template['drmaa_log_dir'] = drmaa_log_dir

    os.makedirs(fna_db_dir, exist_ok=True)
    os.makedirs(bwa_index_dir, exist_ok=True)
    os.makedirs(conda_dir, exist_ok=True)
    os.makedirs(drmaa_log_dir, exist_ok=True)
    if not os.path.exists(config_location):
        os.makedirs(config_location, exist_ok=True)

    config_location = os.path.join(config_location, 'config.yaml')
    with open(config_location, 'w+') as file:
        _ = yaml.dump(config_template, file, sort_keys=False)


    with open(os.path.join(dir_of_this_file, './../config_internal.py'), 'w') as configfile:
        configfile.write('config_loc = "' + config_location + '"')
