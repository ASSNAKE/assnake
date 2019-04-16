from os import path
from os.path import exists
from glob import glob
import json

JAVA = config['java.bin']

FASTQC = config['fastqc.bin']

TRIMMOMATIC = config['trimmomatic']['bin']

def multiple_exts(template, *extensions):
    def wrapped_input(wildcards):
        for ext in extensions:
            p = template.format(**wildcards) + ext
            if exists(p):
                return p
        return p
    return wrapped_input

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
    
    print(params_str)
    
    return params_str



BBDUK = config['bbduk']
rule bbduk_rm_aaa:
    input:
        f="{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz",
        r="{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz",
        rem = 'params/bbduk/aaaaa.fa'
    output:
        f="{prefix}/{df}/reads/{preproc}__bbduk_rmaaa/{sample}/{sample}_R1.fastq.gz",
        r="{prefix}/{df}/reads/{preproc}__bbduk_rmaaa/{sample}/{sample}_R2.fastq.gz",
        stats="{prefix}/{df}/reads/{preproc}__bbduk_rmaaa/{sample}/stats.txt",
    log: "{prefix}/{df}/reads/{preproc}__bbduk_rmaaa/{sample}/{sample}.done"
    shell: ("""{BBDUK} -Xmx1g in1={input.f} in2={input.r} out1={output.f} out2={output.r} ref={input.rem} k=31 hdist=1 stats={output.stats} >{log} 2>&1""")


rule trim_to_len_r2_after_trimmomatic:
    input: 
        r1="datasets/{df}/reads/trimmomatic/{sample}_R1.fastq.gz", 
        r2="datasets/{df}/reads/trimmomatic/{sample}_R2.fastq.gz", 
        u="datasets/{df}/reads/trimmomatic/{sample}_S.fastq.gz"
    output: 
        r1="datasets/{df}/reads/len_trim/{sample}_R1.fastq.gz", 
        r2="datasets/{df}/reads/len_trim/{sample}_R2.fastq.gz", 
        u="datasets/{df}/reads/len_trim/{sample}_S.fastq.gz"
    log: "log/{df}/len_trim/{sample}.trim.log"
    run:
        shell('''cp {input.r1} {output.r1}''')
        shell('''cp {input.u} {output.u}''')
        shell('''{config[python.bin]} bin/scripts/len_trimmer.py --input {input.r2} --out {output.r2} --len 99''')
        
rule fastqc:
    input: "{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz"
    output: 
        zipped="{prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}_fastqc.zip"
    params: 
        out="{prefix}/{df}/reads/{preproc}/{sample}/profile/",
        zip_out="{prefix}/{df}/reads/{preproc}/{sample}/profile/"
    log: "{prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}.log"
    threads: 6
    run:
        shell("{FASTQC} -t {threads} -o {params.out} {input} >{log} 2>&1")
        shell('unzip {output.zipped} -d {params.zip_out}')
        #save_to_db(config['task_id'], rule, str(input), str(output.zipped), 'RUN SUCCESSFUL')

rule count:
    input: 
        r1 = "{fs_prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz",
    output: 
        r1 = "{fs_prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}.count"
    wildcard_constraints:    
        df="[\w\d_-]+"
    run:
        shell('''bin/scripts/count_bp.sh {input.r1} > {output.r1}''')
        #save_to_db(config['task_id'], rule, str(input), str(output.r1), 'RUN SUCCESSFUL')

rule profile:
    input: 
        count = 'datasets/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}.count',
        fastqc = 'datasets/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}_fastqc.zip'
    output:
        profiled = 'datasets/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}_profile.tsv'
    params:
        fastqc_dir='datasets/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}_fastqc/'
    run: 
        import base64
        fastqc_img_dir=params.fastqc_dir+'Images/'
        fastqc_summary=params.fastqc_dir+'summary.txt'
        fastqc_data = params.fastqc_dir + 'fastqc_data.txt'
        result_tsv = ''
        
        count_line = ''
        with open(input.count, 'r') as file:
            for line in file: 
                count_line = line

        reads_str = 'reads\t'+str(int(count_line.split(' ')[0]))
        bp_str = 'bp\t'+str(int(count_line.split(' ')[1]))
        result_tsv += reads_str + '\n' + bp_str + '\n'

        with open(fastqc_summary, 'r') as file:
            for line in file: 
                spl = line.split('\t')
                result_tsv += spl[1].replace(' ', '_').lower() + '\t' + spl[0] + '\n'
                
        
        with open(fastqc_data, 'r') as file:
            for line in file: 
                spl = line.split('\t')
                if (spl[0] == 'File type' or spl[0] == 'Encoding' or 
                    spl[0] == 'Sequences flagged as poor quality' or spl[0] == 'Sequence length'):
                    a = str(spl[0]).replace(' ', '_')
                    a = a.lower()
                    print(a)
                    print(spl[1])
                    line = a + '\t' + spl[1]
                    
                    result_tsv += line
                
        images = ['adapter_content', 'duplication_levels', 'per_base_n_content', 
          'per_base_quality', 'per_base_sequence_content', 'per_sequence_gc_content', 
          'per_sequence_quality', 'sequence_length_distribution']

        for image in images:
            image_loc = fastqc_img_dir + image + '.png'
            result_tsv += image+'_img' + '\t' + image_loc + '\n'
        
        with open(output.profiled, 'w') as file:
            file.write(result_tsv)
            
        if 'task_id' in config.keys():
            save_to_db(config['task_id'], rule, str(input), str(output.profiled), 'RUN SUCCESSFUL')
        
        
        
        
rule fastq_pair:
    input: 
        r1 = "datasets/{df}/reads/{preproc}/{sample}_R1.fastq.gz",
        r2 = "datasets/{df}/reads/{preproc}/{sample}_R2.fastq.gz"
    output: 
        r1 = "datasets/{df}/reads/{preproc}__pair/{sample}_R1.fastq.gz",
        r2 = "datasets/{df}/reads/{preproc}__pair/{sample}_R2.fastq.gz",
    shell:
        '''/data6/bio/TFM/soft/fastq-pair/fastq_pair in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} outsingle={output.s} >{log} 2>&1'''
    
        
rule bbmap_repair:
    input: 
        r1 = "{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz",
        r2 = "{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz"
    output: 
        r1 = "{prefix}/{df}/reads/{preproc}__repair/{sample}/{sample}_R1.fastq.gz",
        r2 = "{prefix}/{df}/reads/{preproc}__repair/{sample}/{sample}_R2.fastq.gz",
        s =  "{prefix}/{df}/reads/{preproc}__repair/{sample}/{sample}_S.fastq.gz"
    log: "{prefix}/{df}/reads/{preproc}__repair/{sample}/{sample}.log"
    shell:
        '''/data6/bio/TFM/soft/bbmap/repair.sh -Xmx16g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} outsingle={output.s} >{log} 2>&1'''

