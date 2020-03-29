import click, os
from assnake.utils import get_config_loc, get_internal_config, pathizer
from assnake.cli.init_config import fill_and_write_config
from assnake.config_internal import config_loc

@click.command(name='start')
@click.option('--assnake_db', '-d', type=click.Path(),
              prompt='Which directory do you like to use for assnake database?',
              help='Desired location of the configiguration file')
@click.option('--fna_db_dir', '-F', type=click.Path(), default='{assnake_db}/fna_db_dir',
              help='This is the place, were all the fasta goes.', required=False)
@click.option('--bwa_index_dir', '-B', type=click.Path(), default='{assnake_db}/bwa_index_dir',
              help='This is the place, were all the bwa index files goes.', required=False)
@click.option('--conda_dir', '-E', type=click.Path(), default='{assnake_db}/conda_dir',
              help='This is the place, were all the conda environments goes.', required=False)
@click.option('--drmaa_log_dir', '-D', type=click.Path(), default='{assnake_db}/drmaa_log_dir',
              help='This is the place, were all the drmaa log files goes.', required=False)
@click.option('--config_location', '-C', type=click.Path(), default='{assnake_db}',
              help='This is the file, were the configure of database will be stored.', required=False)
@click.option('--forced', '-f', is_flag=True,
              help='Force and confirm all in quite mode.')
@click.option('--verbose', '-v', count=True)
# @graph_of_calls('cli_start_init.png')
def init_start(assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, config_location, forced, verbose):
    # TODO Rework assnake install dir
    assnake_db = assnake_db.rstrip('/')

    if config_loc is None or not os.path.isfile(config_loc):
        if not os.path.exists(assnake_db):
            os.mkdir(assnake_db)
        fna_db_dir = fna_db_dir.format(assnake_db=assnake_db)
        bwa_index_dir = bwa_index_dir.format(assnake_db=assnake_db)
        conda_dir = conda_dir.format(assnake_db=assnake_db)
        drmaa_log_dir = drmaa_log_dir.format(assnake_db=assnake_db)
        config_location = config_location.format(assnake_db=assnake_db)
        if verbose > 0:

            def fna_ok():
                return click.confirm(
                    "Proposed folder for fasta files:     " + click.style(fna_db_dir,
                                                                           bold=True,
                                                                           fg='green') + '\n  Is it ok?')

            def bwa_index_ok():
                return click.confirm(
                    "Proposed folder for bwa index files: " + click.style(bwa_index_dir,
                                                                          bold=True,
                                                                          fg='green') + '\n  Is it ok?')

            def conda_ok():
                return click.confirm(
                    "Proposed folder for conda envs:      " + click.style(conda_dir,
                                                                          bold=True,
                                                                          fg='green') + '\n  Is it ok?')

            def log_ok():
                return click.confirm(
                    "Proposed folder for drmaa log files: " + click.style(drmaa_log_dir,
                                                                          bold=True,
                                                                          fg='green') + '\n  Is it ok?')

            def conf_ok():
                return click.confirm(
                    "Proposed path to congig file directory: " + click.style(os.path.join(config_location, 'config.yaml'),
                                                                   bold=True,
                                                                   fg='green') + '\n  Is it ok?')

            while not fna_ok():
                fna_db_dir = click.prompt('Provide absolute path for fasta directory')
            while not bwa_index_ok():
                bwa_index_dir = click.prompt('Provide absolute path for bwa index directory')
            while not conda_ok():
                conda_dir = click.prompt('Provide absolute path for conda environments directory')
            while not log_ok():
                drmaa_log_dir = click.prompt('  Provide absolute path for DRMAA log directory')
            while not conf_ok():
                config_location = click.prompt('  Provide absolute path for configure directory')
            if verbose > 1:
                click.secho(
                    '  Your configuration file will be stored at ' + click.style(os.path.join(config_location, 'config.yaml'), bold=True,
                                                                                 fg='green') + '\n')
        fill_and_write_config(*map(pathizer, [assnake_db, fna_db_dir, bwa_index_dir, conda_dir, drmaa_log_dir, config_location]))
        click.secho('== SUCCESS ==', bg='green', fg='white')
    else:
        if not forced:
            click.echo('We found configuration file at ' + click.style(get_config_loc(), fg='green'))
            if verbose > 0:
                click.secho('Please, keep in mind', err=True)
                click.secho('In (this) case, there is old installation -- only overwriting of ' +
                            click.style('config_internal.ini',
                                        fg='green') + ' in current version is supported.\nDon`t be '
                                                     'upset, cookies are still nice thing and you always can do rm-rf')
            if click.confirm('  Would you like to reset?', abort=True):
                click.echo('As you wish')
        config_internal = get_internal_config()
        dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
        config_internal['GENERAL']['config_loc'] = 'None'
        with open(os.path.join(dir_of_this_file, './../config_internal.ini'), 'w') as configfile:
            config_internal.write(configfile)

@click.command(name = 'current-config')
def current_config():
    """
    Just print your current config location
    """
    click.echo('CURRENT CONFIG LOCATION')
    click.secho(get_config_loc(), bold=True)