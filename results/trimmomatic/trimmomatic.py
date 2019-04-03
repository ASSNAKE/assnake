rule tmtic:
    input:
        first="{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq",
        second="{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq",
        params="params/tmtic/{params}.json"
    output:
        r1="{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}_R1.fastq.gz", 
        r2="{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}_R2.fastq.gz", 
        u ="{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}_S.fastq.gz"
        #                       raw__tmtic_def1__cutadpt
    params:
        u1="{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}_R1_unpaired.fastq", 
        u2="{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}_R2_unpaired.fastq",
    log: "{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}.done"
    threads: 8
    run:
        