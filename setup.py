from setuptools import setup, find_packages
import os
import configparser

# Creating config_internal.ini
dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
config_internal = configparser.ConfigParser()
config_internal['GENERAL'] = {'config_loc': 'None'}
with open(os.path.join(dir_of_this_file, './assnake/config_internal.ini'), 'w+') as configfile:
    config_internal.write(configfile)

setup(name='assnake', 
    version='0.6.0',
    license='MIT',        
    description = 'System for metagenomics data analysis',   
    author = 'Dmitry Fedorov',                   # Type in your name
    author_email = 'fedorov.de@gmail.com',      # Type in your E-Mail
    url = 'https://github.com/Fedorov113/assnake',   # Provide either the link to your github or to your website
    packages=find_packages(),
    install_requires=[
        'numpy', 'Click', 'pyyaml', 'pandas', 
        'tabulate', 'snakemake', 'drmaa', 
        'parse', 'pycallgraph', 'tqdm', 'scipy', 'plotly', 'matplotlib', 'scikit-bio'
    ],
    entry_points='''
        [console_scripts]
        assnake=assnake.cli.assnake_cli:main
    ''',
    )