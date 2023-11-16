import click
import os
from pathlib import Path

from pkg_resources import iter_entry_points
from assnake.core.config import (
    read_internal_config, read_assnake_instance_config,
    update_internal_config, fill_and_write_instance_config
)


def discover_and_deploy_modules():
    discovered_modules = []

    for entry_point in iter_entry_points('assnake.plugins'):
        module_object = entry_point.load()
        discovered_modules.append(module_object)

    for module_object in discovered_modules:
        module_object.deploy_module()


@click.command(name='init')
@click.option('--assnake-db', '-d', type=click.Path(),
              prompt='Directory for Assnake database',
              help='Location to store Assnake database with datasets, tools, and parameters.')
@click.option('--fna-db-dir', '-F', type=click.Path(), default='{assnake_db}/fna_db_dir',
              help='Directory for FASTA files.', required=False)
@click.option('--bwa-index-dir', '-B', type=click.Path(), default='{assnake_db}/bwa_index_dir',
              help='Directory for BWA index files.', required=False)
@click.option('--conda-dir', '-E', type=click.Path(), default='{assnake_db}/conda_dir',
              help='Directory for Conda environments.', required=False)
@click.option('--drmaa-log-dir', '-D', type=click.Path(), default='{assnake_db}/drmaa_log_dir',
              help='Directory for DRMAA log files.', required=False)
@click.option('--start-over', is_flag=True, default=False,
              help='Reconfigure Assnake from scratch.')
@click.option('--forced', '-f', is_flag=True,
              help='Force all confirmations in quiet mode.')
@click.option('--verbose', '-v', count=True)
def init_config(assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, start_over, forced, verbose):
    """Initialize Assnake installation."""
    internal_config = read_internal_config()
    assnake_db = str(Path(assnake_db).expanduser().resolve())
    instance_config_loc = os.path.join(assnake_db, 'config.yaml')

    if start_over or read_assnake_instance_config() is None:
        if not os.path.isfile(instance_config_loc):
            os.makedirs(assnake_db, exist_ok=True)
            fna_db_dir = fna_db_dir.format(assnake_db=assnake_db)
            bwa_index_dir = bwa_index_dir.format(assnake_db=assnake_db)
            conda_dir = conda_dir.format(assnake_db=assnake_db)
            drmaa_log_dir = drmaa_log_dir.format(assnake_db=assnake_db)
            
            # Prompt user for inputs if not forced and verbose
            if not forced and verbose:
                fna_db_dir = click.prompt('Path for FASTA database', default=fna_db_dir)
                bwa_index_dir = click.prompt('Path for index files', default=bwa_index_dir)
                conda_dir = click.prompt('Path for Conda environments', default=conda_dir)
                drmaa_log_dir = click.prompt('Path for log files', default=drmaa_log_dir)
            
            instance_config_loc = fill_and_write_instance_config(
                assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, instance_config_loc
            )

            if instance_config_loc:
                update_internal_config({'instance_config_loc': instance_config_loc})
                # Discover and deploy modules during initialization
                discover_and_deploy_modules()
                click.secho('Assnake database initialized at: ' + assnake_db, fg='green')
        else:
            click.secho('Existing Assnake instance found at: ' + assnake_db, fg='yellow')
            internal_config = update_internal_config({'instance_config_loc': instance_config_loc})
    else:
        click.secho('Assnake already configured. Database at: ' + assnake_db, fg='yellow')

@click.command(name='show-internal')
def show_internal_config():
    """Show the current internal configuration."""
    click.echo(read_internal_config())

@click.command(name='show-instance')
def show_instance_config():
    """Show the current instance configuration."""
    click.echo(read_assnake_instance_config())
