from snakemake.shell import shell
import json
import yaml

def bwa_params(params_loc):
    params_dict = {}
    with open(params_loc, 'r') as stream:
        try:
            params_dict = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
        
    
    
    
    return params_str

param_str = tmtic_params(snakemake.input.params)
        
shell('''trimmomatic PE -phred33 \
                 -threads {snakemake.threads} \
                 {snakemake.input.first} {snakemake.input.second} \
                 {snakemake.output.r1} {snakemake.params.u1} \
                 {snakemake.output.r2} {snakemake.params.u2} \
                 {param_str} \
         >{snakemake.log} 2>&1 && \
         cat {snakemake.params.u1} {snakemake.params.u2} | gzip > {snakemake.output.u} 2>>{snakemake.log} && \
         rm {snakemake.params.u1} {snakemake.params.u2} 2>>{snakemake.log}''')

if 'task_id' in snakemake.config.keys():
    save_to_db(config['task_id'], 'tmtic', str(input), str(log), 'RUN SUCCESSFUL')