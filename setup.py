from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

import yaml
import fastentrypoints
import os
import configparser
from pathlib import Path
# Creating config_internal.ini
# dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
# config_internal = configparser.ConfigParser()
# config_internal['GENERAL'] = {'config_loc': 'None'}
# with open(os.path.join(dir_of_this_file, './assnake/config_internal.ini'), 'w+') as configfile:
#     config_internal.write(configfile)

def write_internal_config():
    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')
    internal_config_dir = os.path.join(str(Path.home()), '.config/assnake/')

    os.makedirs(internal_config_dir, exist_ok=True)

    if not os.path.isfile(internal_config_loc):
        with open(internal_config_loc, 'w+') as file:
            _ = yaml.dump({'instance_config_loc': 'not_set'}, file, sort_keys=False)

class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        write_internal_config()
        develop.run(self)

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        write_internal_config()
        install.run(self)

setup(name='assnake', 
    version='0.9.0.1',
    include_package_data=True,
    license='MIT',         
    description = 'System for metagenomics data analysis',   
    author = 'Dmitry Fedorov',                  
    author_email = 'fedorov.de@gmail.com',      
    url = 'https://github.com/ASSNAKE/assnake',   
    download_url = 'https://github.com/ASSNAKE/assnake/archive/v0.8.8.tar.gz',    # I explain this later on
    keywords = ['ILLUMINA', 'NGS', 'METAGENOMIC', 'DATA'], 
    packages=find_packages(),
    install_requires=[
        'numpy', 'Click', 'pyyaml>=5', 'pandas', 
        'tabulate', 'snakemake', 'drmaa', 
        'parse', 'pycallgraph', 'tqdm', 'scipy', 'plotly', 'matplotlib'
    ],
    entry_points='''
        [console_scripts]
        assnake=assnake.cli.assnake_cli:main
    ''',

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Science/Research',      # Define that your audience are developers
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',   # Again, pick a license
        'Programming Language :: Python :: 3.6',
        ],
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
        }
    )