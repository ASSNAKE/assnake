from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

import fastentrypoints
import os
import configparser
from pathlib import Path

def write_internal_config():
    internal_config_loc = os.path.join(str(Path.home()), '.config/assnake/internal_config.yaml')
    internal_config_dir = os.path.join(str(Path.home()), '.config/assnake/')

    os.makedirs(internal_config_dir, exist_ok=True)

    if not os.path.isfile(internal_config_loc):
        with open(internal_config_loc, 'w+') as file:
            file.write('instance_config_loc: not_set')

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
    version='0.9.0.2',
    include_package_data=True,
    license='MIT',         
    description = 'System for metagenomics data analysis',   
    author = 'Dmitry Fedorov',                  
    author_email = 'fedorov.de@gmail.com',      
    url = 'https://github.com/ASSNAKE/assnake',   
    # download_url = 'https://github.com/ASSNAKE/assnake/archive/v0.8.8.tar.gz',    # I explain this later on
    keywords = ['ILLUMINA', 'NGS', 'METAGENOMIC', 'DATA'], 
    packages=find_packages(),
    install_requires=[
        'numpy', 'Click', 'pyyaml', 'pandas', 
        'tabulate', 'snakemake', 'tqdm', 'scipy', 'matplotlib', 'plotly', 'parse'
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
        'Programming Language :: Python :: 3.10',
        ],
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
        }
    )
