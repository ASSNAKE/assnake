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


curr_dir = os.path.abspath(os.path.dirname(__file__))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, default='',
                        help='Snakefile with defenitions of your request')
    parser.add_argument('-n', '--dry-run', action='store_true')
    args = parser.parse_args()


    cwd = os.getcwd()
    print(args.s)
    print(curr_dir)

    snakefile = args.s
    status = snakemake.snakemake(snakefile, 
                                config = dict(assnake_install_dir= curr_dir),    
                                targets=['test_count'], 
                                printshellcmds=True,
                                dryrun=args.dry_run, 
                                configfile=os.path.join(curr_dir, 'config.yml'),
                                use_conda = True
                                )

    if 1: # translate "success" into shell exit code of 0
       return 0
    return 1



if __name__ == '__main__':
    main()