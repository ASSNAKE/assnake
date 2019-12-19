from snakemake.shell import shell
import yaml
import os

def get_params_str(params_loc):
    params_str = '{truncLen_f} {truncLen_r} {trimLeft_l} {trimLeft_r} {maxEE_f} {maxEE_r} {truncQ} {maxN}'
    with open(snakemake.input.params, 'r') as stream:
        try:
            par = yaml.safe_load(stream)
            params_str = params_str.format(**par)
        except yaml.YAMLError as exc:
            print(exc)
    
    return params_str

params_str = get_params_str(snakemake.input.params)
filter_trim_script = os.path.join(snakemake.config['assnake_install_dir'], 'modules/dada2/scripts/filter_trim.R')
shell('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
                    Rscript {filter_trim_script} '{snakemake.input.r1}' '{snakemake.input.r2}' '{snakemake.output.r1}' '{snakemake.output.r2}' {params_str} >{snakemake.log} 2>&1''')

# if 'task_id' in snakemake.config.keys():
#     save_to_db(config['task_id'], 'tmtic', str(input), str(log), 'RUN SUCCESSFUL')