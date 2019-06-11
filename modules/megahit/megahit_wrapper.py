from snakemake.shell import shell
import json
import shutils

def megahit_params(params_loc):
    params_str = ''
    params_dict = {}
    
    return params_str

# param_str = megahit_params(snakemake.input.params)
        
print(snakemake.input.R)
print(snakemake.input.F)

if os.path.exists(snakemake.params.out_folder) and os.path.isdir(snakemake.params.out_folder):
    shutil.rmtree(snakemake.params.out_folder)

shell('''megahit -1 {reads1} -2 {reads2} --min-contig-len 850 -o {snakemake.params.out_folder} -t {snakemake.threads}  >{snakemake.log} 2>&1''')

if 'task_id' in snakemake.config.keys():
    save_to_db(config['task_id'], 'tmtic', str(input), str(log), 'RUN SUCCESSFUL')