from snakemake.shell import shell
import json
import shutil
import os

def megahit_params(params_loc):
    params_dict = {}
    params_str = ''
    with open(params_loc, 'r') as stream:
        params_dict = json.load(stream)

    bao = params_dict['basic_assembly_options']
    bao_wc = '--min-count {min_count} --k-list {k_list} '
    params_str += bao_wc.format(min_count = str(bao['min_count']), k_list=','.join((str(x) for x in bao['k_list'])))

    aao = params_dict['advanced_assembly_options']
    aao_wc = '--bubble-level {bubble_level} --prune-level {prune_level} --prune-depth {prune_depth} --low-local-ratio {low_local_ratio} --max-tip-len {max_tip_len} --merge-level {merge_level_l},{merge_level_s} '
    if aao['no_mercy']:
        params_str += '--no-mercy '
    if aao['no_local']:
        params_str += '--no-local '
    if aao['kmin_1pass']:
        params_str += '--kmin-1pass '
    params_str += aao_wc.format(bubble_level = aao['bubble_level'],
                            prune_level = aao['prune_level'],
                            prune_depth = aao['prune_depth'],
                            low_local_ratio = aao['low_local_ratio'],
                            max_tip_len = aao['max_tip_len'],
                            merge_level_l = aao['merge_level']['merge_level_l'],
                            merge_level_s = aao['merge_level']['merge_level_s'])

    params_str += '--min-contig-len {}'.format(params_dict['output_options']['min_contig_len'])
    
    return params_str

params_str = megahit_params(snakemake.input.params)
reads1 = ','.join(snakemake.input.F)
reads2 = ','.join(snakemake.input.R)

shell('''rm -rf {snakemake.params.out_folder}; megahit -1 {reads1} -2 {reads2} {params_str} -o {snakemake.params.out_folder} -t {snakemake.threads} > {snakemake.log} 2>&1''')
fc_loc = snakemake.params.out_folder+'final.contigs.fa'
shell('cp {fc_loc} {snakemake.output.out_fa}; cp {snakemake.input.params} {snakemake.output.params}')

if 'task_id' in snakemake.config.keys():
    save_to_db(config['task_id'], 'tmtic', str(input), str(log), 'RUN SUCCESSFUL')