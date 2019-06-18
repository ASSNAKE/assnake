from snakemake.shell import shell
import yaml
import os

def get_params_str(params_loc):
    params_str = "'{randomize}' {MAX_CONSIST}"
    with open(snakemake.input.params, 'r') as stream:
        try:
            par = yaml.safe_load(stream)
            params_str = params_str.format(**par)
        except yaml.YAMLError as exc:
            print(exc)
    
    return params_str

params_str = get_params_str(snakemake.input.params)
learn_errors_script = os.path.join(snakemake.config['assnake_install_dir'], 'modules/dada2/scripts/learn_errors.R')
shell('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript  {learn_errors_script} '{snakemake.input.samples_list}' \
            '{snakemake.output.err}' '{snakemake.wildcards.strand}' {params_str} > {snakemake.log} 2>&1''')
