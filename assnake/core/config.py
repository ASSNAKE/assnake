import os
from pathlib import Path
import yaml
from assnake.utils.general import read_yaml
import click

# Constants for file paths
INTERNAL_CONFIG_FILE = Path.home() / '.config/assnake/internal_config.yaml'
CONFIG_TEMPLATE_FILE = Path(__file__).parent / '../snake/config_template.yml'


def load_wc_config():
    """
    Load the wildcard configuration for snakemake workflows.
    """
    return read_yaml(Path(__file__).parent / '../snake/wc_config.yaml')


def read_internal_config():
    """
    Read the internal configuration from the user's home directory.

    Returns:
        A dictionary representing the internal configuration.
    """
    if not INTERNAL_CONFIG_FILE.is_file():
        os.makedirs(INTERNAL_CONFIG_FILE.parent, exist_ok=True)
        with open(INTERNAL_CONFIG_FILE, 'w') as file:
            yaml.dump({'instance_config_loc': 'not_set'}, file, sort_keys=False)
        return {'instance_config_loc': 'not_set'}

    return read_yaml(INTERNAL_CONFIG_FILE)


def read_assnake_instance_config():
    """
    Read the instance configuration for the assnake system.

    Returns:
        A dictionary of the instance configuration if the file exists, None otherwise.
    """
    internal_config = read_internal_config()
    instance_config_loc = internal_config.get('instance_config_loc')

    return read_yaml(instance_config_loc) if instance_config_loc and Path(instance_config_loc).is_file() else None


def update_config(config_loc, config_dict):
    """
    Generic function to update a YAML configuration file.

    Args:
        config_loc: Location of the configuration file.
        config_dict: Dictionary with configuration to be updated.
    """
    with open(config_loc, 'w') as file:
        yaml.dump(config_dict, file, sort_keys=False)


def update_internal_config(update_dict):
    """
    Update the internal configuration with the provided dictionary.
    """
    internal_config = read_internal_config()
    internal_config.update(update_dict)
    update_config(INTERNAL_CONFIG_FILE, internal_config)
    return internal_config


def check_if_assnake_is_initialized():
    """
    Check if the assnake system has been initialized by reading the instance config.
    If not initialized, prompts the user to initialize.
    """
    if read_assnake_instance_config() is None:
        click.secho("Assnake is not initialized. Run `assnake config init` to set up.", fg='red', bold=True)
        exit()


def fill_and_write_instance_config(assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, instance_config_loc):
    """
    Fill and write the instance configuration based on provided locations.

    Args:
        assnake_db: Database directory for assnake.
        fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir: Directories for respective data types.
        instance_config_loc: Location to save the instance configuration file.

    Returns:
        The location of the instance configuration file if successful, None otherwise.
    """
    config_template = read_yaml(CONFIG_TEMPLATE_FILE)
    config_template.update({
        'assnake_db': assnake_db,
        'assnake_install_dir': Path(__file__).parent.parent,  # updated to use relative path
        'fna_db_dir': fna_db_dir,
        'bwa_index_dir': bwa_index_dir,
        'conda_dir': conda_dir,
        'drmaa_log_dir': drmaa_log_dir
    })

    for dir_path in [fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir]:
        os.makedirs(dir_path, exist_ok=True)

    update_config(instance_config_loc, config_template)
    return instance_config_loc