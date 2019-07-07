from setuptools import setup, find_packages

setup(name='assnake', 
    version='0.0.3', 
    packages=find_packages(),
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        assnake=assnake:cli
    ''',
    )