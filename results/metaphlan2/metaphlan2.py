rule metaphlan2:
    input:
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
    output:
        o = '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/{sample}.mp2'
    params:
        b =  '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/{sample}.b2',
        MPA_PKL = config['MetaPhlAn2']['mpa_v29'],
        BOWTIE2DB = config['MetaPhlAn2']['mpa_v29_bt'],
        INDEX = 'v29_CHOCOPhlAn_201901',
        task_id = config['task_id'] if 'task_id' in config.keys() else None,
    log: '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/log.txt'
    benchmark: '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/time.txt'
    threads: 12
    conda: 'env_2.7.8.yaml'
    wrapper: "file:///data6/bio/TFM/pipeline/assnake/results/metaphlan2/wrapper.py"
        
