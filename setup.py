from setuptools import setup, find_packages

setup(name='assnake', 
    version='0.3.0',
    license='MIT',        
    description = 'System for metagenomics data analysis',   
    author = 'Dmitry Fedorov',                   # Type in your name
    author_email = 'fedorov.de@gmail.com',      # Type in your E-Mail
    url = 'https://github.com/Fedorov113/assnake',   # Provide either the link to your github or to your website
    packages=find_packages(),
    install_requires=[
        'Click', 'pyyaml', 'pandas', 'tabulate', 'snakemake', 'drmaa'
    ],
    entry_points='''
        [console_scripts]
        assnake=assnake.assnake_cli:main
    ''',
    )