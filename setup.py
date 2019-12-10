from setuptools import setup, find_packages

setup(name='assnake', 
    version='0.0.9', 
    packages=find_packages(),
    install_requires=[
        'Click', 'pyyaml', 'pandas', 'tabulate', 'snakemake', 'drmaa'
    ],
    entry_points='''
        [console_scripts]
        assnake=assnake.assnake_cli:cli
    ''',
    )