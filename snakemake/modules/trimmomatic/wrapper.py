from snakemake.shell import shell
import json

def tmtic_params(params_loc):
    params_str = ''
    params_dict = {}
    with open(params_loc, 'r') as f:
        params_dict = json.loads(f.read())
        
    clip = params_dict.pop('ILLUMINACLIP')
    if(len(clip.keys()) > 0):
        slw_str = 'ILLUMINACLIP:{fastaWithAdaptersEtc}:{seedMismatches}:{palindromeClipThreshold}:{simpleClipThreshold} '
        params_str += slw_str.format(**clip)
    
    slw = params_dict.pop('SLIDINGWINDOW')
    if(len(slw.keys()) > 0):
        slw_str = 'SLIDINGWINDOW:{windowSize}:{requiredQuality} '
        params_str += slw_str.format(**slw)
        
    lead = params_dict.pop('LEADING')
    if(len(lead.keys()) > 0):
        slw_str = 'LEADING:{quality} '
        params_str += slw_str.format(**lead)
    
    trail = params_dict.pop('TRAILING')
    if(len(trail.keys()) > 0):
        slw_str = 'TRAILING:{quality} '
        params_str += slw_str.format(**trail)
        
    minlen = params_dict.pop('MINLEN')
    if(len(minlen.keys()) > 0):
        slw_str = 'MINLEN:{length} '
        params_str += slw_str.format(**minlen)
        
    hcrop = params_dict.pop('HEADCROP')
    if(len(hcrop.keys()) > 0):
        slw_str = 'HEADCROP:{length} '
        params_str += slw_str.format(**hcrop)
    
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