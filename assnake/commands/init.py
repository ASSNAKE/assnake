import click, sys, os, glob, yaml, configparser

from assnake.utils import read_yaml

@click.command(name='start')
@click.option('--config_location','-c', type=click.Path(), 
    # prompt='Desired location of the config file', 
    help='Desired location of the configiguration file' )
@click.option('--assnake_db','-c', type=click.Path(), 
    # prompt='Directory for assnake database', 
    help='Assnake database', required=False)
@click.option('--fna_db_dir','-c', type=click.Path(), 
    # prompt='Directory for your fasta files', 
    help='This is the place, were all the fasta goes.' )
def init_start(config_location, assnake_db, fna_db_dir):
    overwrite = False

    dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
    config_internal = configparser.ConfigParser()
    config_internal.read(os.path.join(dir_of_this_file, './../config_internal.ini'))
    config_loc = config_internal['GENERAL']['config_loc']

    config_template_loc = os.path.join(dir_of_this_file, '../../snake/config_template.yml')
    config_template = read_yaml(config_template_loc)
    # TODO Rework assnake install dir
    assnake_install_dir = (dir_of_this_file[0:-7]+'snake')
    if config_loc == 'None':
        click.echo(click.style('Which directory do you want to use for assnake database', fg='green'))
        assnake_db = click.prompt('  Provide absolute path to the directory')
        assnake_db = assnake_db.replace(' ', '')
        if (os.path.exists(assnake_db)):
            use_folder_as_assnake_db = click.confirm("  Folder exists! Do you want to use it?")
            if use_folder_as_assnake_db:
                config_location = os.path.join(assnake_db, 'config.yaml')
                if not os.path.isfile(config_location) or True:
                    config_template['assnake_db'] = assnake_db
                    config_template['assnake_install_dir'] = assnake_install_dir

                    click.secho('  Your configuration file will be stored at ' +  click.style(config_location, bold = True, bg='blue') +'\n')
                    
                    fna_db_dir = os.path.join(assnake_db, 'fna_db_dir')
                    bwa_index_dir = os.path.join(assnake_db, 'bwa_index_dir')
                    conda_dir = os.path.join(assnake_db, 'conda_dir')
                    drmaa_log_dir = os.path.join(assnake_db, 'drmaa_log_dir')

                    fna_ok       = click.confirm("  Proposed folder for fasta files:     " + click.style(fna_db_dir, bold = True, bg='blue') +' Is it ok?\n')
                    if not fna_ok:
                        fna_db_dir = click.prompt('  Provide absolute path for fasta directory')

                    bwa_index_ok = click.confirm("  Proposed folder for bwa index files: " + click.style(bwa_index_dir, bold = True, bg='blue') +' Is it ok?\n')
                    if not bwa_index_ok:
                        bwa_index_dir = click.prompt('  Provide absolute path for bwa index directory')

                    conda_ok     = click.confirm("  Proposed folder for conda envs:      " + click.style(conda_dir, bold = True, bg='blue') +' Is it ok?\n')
                    if not conda_ok:
                        conda_dir = click.prompt('  Provide absolute path for conda environments directory')

                    log_ok       = click.confirm("  Proposed folder for drmaa log files: " + click.style(drmaa_log_dir, bold = True, bg='blue') +' Is it ok?\n')
                    if not log_ok:
                        drmaa_log_dir = click.prompt('  Provide absolute path for DRMAA log directory')

                    config_template['fna_db_dir'] = fna_db_dir
                    config_template['bwa_index_dir'] = bwa_index_dir
                    config_template['conda_dir'] = conda_dir
                    config_template['drmaa_log_dir'] = drmaa_log_dir

                    os.makedirs(fna_db_dir, exist_ok = True)
                    os.makedirs(bwa_index_dir, exist_ok = True)
                    os.makedirs(conda_dir, exist_ok = True)
                    os.makedirs(drmaa_log_dir, exist_ok = True)

                    with open(config_location, 'w') as file:
                        documents = yaml.dump(config_template, file, sort_keys=False)

                    config_internal['GENERAL']['config_loc'] = config_location
                    with open(os.path.join(dir_of_this_file, './../config_internal.ini'), 'w') as configfile:
                        config_internal.write(configfile)
    else:
        click.echo('We found configuration file at ' + click.style(config_loc, bg='blue'))  
    #     click.echo(click.style('Where to store your configuration file?', fg='green'))
    #     config_location = click.prompt('  Provide absolute location of the configuration file')
    #     config_location = config_location.replace(' ', '')
    #     if (os.path.exists(config_location)):
    #         click.echo("This file exists!")
    #     else:
    #         config_internal['GENERAL']['config_loc'] = config_location
    #         try:
    #             os.makedirs(os.path.dirname(config_location), exist_ok=True)
    #             shutil.copyfile(os.path.join(dir_of_this_file, '../snake/config_template.yml'), config_location)
    #             click.secho('Config file initialized in ' + config_location, fg='green')
    #         except:
    #             click.secho('ERROR CREATING '+os.path.dirname(config_location) , fg='red')
            
    #         with open(os.path.join(dir_of_this_file, './config_internal.ini'), 'w') as configfile:
    #            config_internal.write(configfile)
    # else:
    #     click.secho('Configuration file: ' + config_internal['GENERAL']['config_loc'])
    #     config_loc = config_internal['GENERAL']['config_loc']
    #     config = read_yaml(config_loc)
    #     pprint.pprint(config)
        # print(config_loc)



    # config_loc = os.path.join(dir_of_this_file, '../snake/config.yml')

    

    # # if not os.path.isdir(config_location):
    # #     print(config_location)
    # #     os.makedirs(config_location, exist_ok=True)

    # # with open(os.path.join(config_location, 'config.yml'), 'w') as config_file:
    # #     written_config = yaml.dump(config, config_file)

    # assnake_db = click.prompt('Please, enter desired location of the ASSNAKE database')
    # config['assnake_db'] = assnake_db
    # config['assnake_install_dir'] = os.path.join(os.path.dirname(dir_of_this_file), 'snakemake')
    
    # with open(config_loc, 'w') as config_file:
    #     written_config = yaml.dump(config, config_file)

    # # shutil.copyfile(config_loc, os.path.join(config_location, 'config.yml'))