import yaml

rule dada2_filater_and_trim:
    input: 
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        params = os.path.join(config['assnake_db'], 'params/dada2/filter_and_trim/{params}.yaml')
    output:
        r1 = '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}_R2.fastq.gz'
    log: '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}.log'
    conda: 'dada2.yaml'
    wrapper: "file:///data6/bio/TFM/pipeline/assnake/results/dada2/filter_trim_wrapper.py"


rule dada2_derep_dada_merge:
    input: 
        r1     = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2     = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        errF   = '/data6/bio/TFM/dada2/{run}/errF.rds',
        errR   = '/data6/bio/TFM/dada2/{run}/errR.rds'
    output:
        merged = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}__{run}.merged.rds',
    log: '{prefix}/{df}/reads/{preproc}/{sample}/{sample}__{run}.log'
    conda: 'dada2.yaml'
    shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript  /data6/bio/TFM/pipeline/assnake/results/dada2/scripts/derep_dada_merge.R '{input.r1}' '{input.r2}' '{input.errF}' '{input.errR}' '{output.merged}' 8 >{log} 2>&1''')