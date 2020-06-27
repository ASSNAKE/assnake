import click, os

from assnake.core.config import read_internal_config, read_assnake_instance_config, update_internal_config, fill_and_write_instance_config

from assnake.core.result import get_all_results_as_list

from pathlib import Path

@click.command(name='init')
@click.option('--assnake-db', '-d', type=click.Path(),
              prompt='Which directory do you like to use for assnake database?',
              help='Assnake database stores information about your datasets, installed tools and parameters.')

@click.option('--fna_db_dir', '-F', type=click.Path(), default='{assnake_db}/fna_db_dir',
              help='This is the place, were all the fasta goes.', required=False)
@click.option('--bwa_index_dir', '-B', type=click.Path(), default='{assnake_db}/bwa_index_dir',
              help='This is the place, were all the bwa index files goes.', required=False)
              
@click.option('--conda_dir', '-E', type=click.Path(), default='{assnake_db}/conda_dir',
              help='This is the place, were all the conda environments goes.', required=False)
@click.option('--drmaa_log_dir', '-D', type=click.Path(), default='{assnake_db}/drmaa_log_dir',
              help='This is the place, were all the drmaa log files goes.', required=False)

@click.option('--start-over', is_flag=True, default = False,
              help='Configure assnake from scratch')

@click.option('--forced', '-f', is_flag=True,
              help='Force and confirm all in quite mode.')
@click.option('--verbose', '-v', count=True)

def init_config(assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir,start_over, forced, verbose):
    "Initialize your Assnake installation"
    internal_config = read_internal_config()

    assnake_db = str(Path(assnake_db).expanduser().resolve())

    # Check if assnake_db in internal config is valid by checking config.yaml
    if read_assnake_instance_config() is None or start_over: # No valid config in internal config
        # Check if instance config file exists in specified assnake_db folder
        if not os.path.isfile(os.path.join(assnake_db, 'config.yaml')): # No instance config file
            # Create dir if not created
            if not os.path.exists(assnake_db):
                os.mkdir(assnake_db)
            # Generate proposed locations
            fna_db_dir          = fna_db_dir.format(assnake_db=assnake_db)
            bwa_index_dir       = bwa_index_dir.format(assnake_db=assnake_db)
            conda_dir           = conda_dir.format(assnake_db=assnake_db)
            drmaa_log_dir       = drmaa_log_dir.format(assnake_db=assnake_db)
            instance_config_loc = os.path.join(assnake_db, 'config.yaml')

            if verbose > 0:

                def fna_ok():
                    return click.confirm("Proposed folder for fasta files: " + click.style(fna_db_dir, bold=True, fg='green') + '\n  Is it ok?')

                def bwa_index_ok():
                    return click.confirm("Proposed folder for index files: " + click.style(bwa_index_dir, bold=True, fg='green') + '\n  Is it ok?')

                def conda_ok():
                    return click.confirm("Proposed folder for conda envs:  " + click.style(conda_dir, bold=True, fg='green') + '\n  Is it ok?')

                def log_ok():
                    return click.confirm("Proposed folder for log files:   " + click.style(drmaa_log_dir, bold=True, fg='green') + '\n  Is it ok?')

                while not fna_ok():
                    fna_db_dir = click.prompt('Provide absolute path for fasta database directory')
                while not bwa_index_ok():
                    bwa_index_dir = click.prompt('Provide absolute path for index directory')
                while not conda_ok():
                    conda_dir = click.prompt('Provide absolute path for conda environments directory')
                while not log_ok():
                    drmaa_log_dir = click.prompt('  Provide absolute path for log directory')

                if verbose > 1:
                    click.secho(
                        '  Your configuration file will be stored at ' + click.style(os.path.join(config_location, 'config.yaml'), bold=True,
                                                                                    fg='green') + '\n')
            instance_config_loc = fill_and_write_instance_config(*map(os.path.abspath, [assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, os.path.join(assnake_db, 'config.yaml')]))

            if instance_config_loc is not None:
                internal_config = update_internal_config({ 'instance_config_loc': instance_config_loc })
                if internal_config is not None:
                    click.secho('We created assnake database instance at ' + assnake_db, fg='green')
                    click.secho('== SUCCESS ==', bg='green', fg='white')
        else:# We found existing config file
            click.secho('We found existing assnake instance configuration at ' + assnake_db, fg='yellow')
            click.secho('We will set it as assnake instance configuration with no changes', fg='yellow')
            internal_config = update_internal_config({
                    'instance_config_loc': os.path.join(assnake_db, 'config.yaml'),
                    })
            if internal_config is not None:
                click.secho('== SUCCESS ==', bg='green', fg='white')

    else:# We are already pointed to a valid assnake instance configuration
        click.secho('We are already properly configured, assnake database is here: ' + assnake_db, fg='yellow')
        click.secho('If you want to start over, add flag --start-over', fg='yellow')

@click.command(name = 'show-internal')
def show_internal_config():
    """
    Show your current internal configuration
    """
    # click.echo('CURRENT CONFIG LOCATION')
    print(read_internal_config())

@click.command(name = 'show-instance')
def show_instance_config():
    """
    Show your current instance configuration
    """
    # click.echo('CURRENT CONFIG LOCATION')
    print(read_assnake_instance_config())

