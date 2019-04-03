
"""
Command line interface driver for snakemake workflows
"""
import argparse
import os.path
import snakemake
import sys
import pprint
import json



# thisdir = os.path.abspath(os.path.dirname(__file__))
# parentdir = os.path.join(thisdir,'..')
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(description='bananas: run snakemake workflows', usage='''bananas <workflow> <parameters> [<target>]
bananas: run snakemake workflows, using the given workflow name & parameters file.
''')

    parser.add_argument('workflowfile')
    parser.add_argument('paramsfile')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-f', '--force', action='store_true')
    args = parser.parse_args(sysargs)

    # first, find the Snakefile
    snakefile_this      = os.path.join(thisdir,"Snakefile")
    snakefile_parent    = os.path.join(parentdir,"Snakefile")
    if os.path.exists(snakefile_this):
        snakefile = snakefile_this
    elif os.path.exists(snakefile_parent):
        snakefile = snakefile_parent
    else:
        msg = 'Error: cannot find Snakefile at any of the following locations:\n'
        msg += '{}\n'.format(snakefile_this)
        msg += '{}\n'.format(snakefile_parent)
        sys.stderr.write(msg)
        sys.exit(-1)

    # next, find the workflow config file
    workflowfile = None
    w1 = os.path.join(cwd,args.workflowfile)
    w2 = os.path.join(cwd,args.workflowfile+'.json')
    if os.path.exists(w1) and not os.path.isdir(w1):
        workflowfile = w1
    elif os.path.exists(w2) and not os.path.isdir(w2):
        workflowfile = w2

    if not workflowfile:
        msg = 'Error: cannot find workflowfile {} or {} '.format(w1,w2)
        msg += 'in directory {}\n'.format(cwd)
        sys.stderr.write(msg)
        sys.exit(-1)

    # next, find the workflow params file
    paramsfile = None
    p1 = os.path.join(cwd,args.paramsfile)
    p2 = os.path.join(cwd,args.paramsfile+'.json')
    if os.path.exists(p1) and not os.path.isdir(p1):
        paramsfile = p1
    elif os.path.exists(p2) and not os.path.isdir(p2):
        paramsfile = p2

    if not paramsfile:
        msg = 'Error: cannot find paramsfile {} or {} '.format(p1,p2)
        msg += 'in directory {}\n'.format(cwd)
        sys.stderr.write(msg)
        sys.exit(-1)

    with open(workflowfile, 'rt') as fp:
        workflow_info = json.load(fp)

    target = workflow_info['workflow_target']
    config = dict()

    print('--------')
    print('details!')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfig: {}'.format(workflowfile))
    print('\tparams: {}'.format(paramsfile))
    print('\ttarget: {}'.format(target))
    print('--------')

    # run bananas!!
    status = snakemake.snakemake(snakefile, configfile=paramsfile,
                                 targets=[target], printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force,
                                 config=config)

    if status: # translate "success" into shell exit code of 0
        return 0
    return 1


if __name__ == '__main__':
    main()