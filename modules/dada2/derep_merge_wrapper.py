from snakemake.shell import shell
import yaml
import os

def get_params_str(params_loc):
    params_str = '{minOverlap}'
    with open(snakemake.input.params, 'r') as stream:
        try:
            par = yaml.safe_load(stream)
            params_str = params_str.format(**par)
        except yaml.YAMLError as exc:
            print(exc)
    return params_str

params_str = get_params_str(snakemake.input.params)
derep_dada_merge_script = os.path.join(snakemake.config['assnake_install_dir'], 'modules/dada2/scripts/derep_dada_merge.R')
shell('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
                    Rscript  {derep_dada_merge_script} '{snakemake.input.r1}' '{snakemake.input.r2}'\
                         '{snakemake.input.errF}' '{snakemake.input.errR}' '{snakemake.output.merged}'\
                              {params_str} '{snakemake.output.stats}' >{snakemake.log} 2>&1''')
