"""
Command line interface driver for ASSNAKE
"""

import argparse
import os.path
import os
import snakemake
import sys
import pprint
import json
import yaml

curr_dir = os.path.abspath(os.path.dirname(__file__))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, default='',
                        help='Snakefile with defenitions of your request')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument("target",
                        nargs="*",
                        default=None,
                        help="Targets to build. May be rules or files.")
    args = parser.parse_args()


    cwd = os.getcwd()
    print(args.s)
    print(curr_dir)

    configfile=os.path.join(curr_dir, 'config.yml')
    conda_prefix = ''
    with open(configfile) as f:
        config = yaml.load(f)
        conda_prefix = config['conda_prefix']

    print(args.target)
    snakefile = args.s
    status = snakemake.snakemake(snakefile, 
                                config = dict(assnake_install_dir= curr_dir),    
                                targets=args.target, 
                                printshellcmds=True,
                                dryrun=args.dry_run, 
                                configfile=os.path.join(curr_dir, 'config.yml'),
                                use_conda = True,
                                conda_prefix = conda_prefix
                                )

    if 1: # translate "success" into shell exit code of 0
       return 0
    return 1



if __name__ == '__main__':
    main()