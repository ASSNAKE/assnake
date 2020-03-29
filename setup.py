from setuptools import setup, find_packages
import fastentrypoints
import os
import configparser

# Creating config_internal.ini
# dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
# config_internal = configparser.ConfigParser()
# config_internal['GENERAL'] = {'config_loc': 'None'}
# with open(os.path.join(dir_of_this_file, './assnake/config_internal.ini'), 'w+') as configfile:
#     config_internal.write(configfile)

setup(name='assnake', 
    version='0.8.9.8',
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
    )