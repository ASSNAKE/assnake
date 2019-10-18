rule kraken2_silva:
    input: 
        r1  = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2  = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
    output: 
        f   = '{prefix}/{df}/taxa/kraken2/{preproc}/{sample}_1.fq',
        s   = '{prefix}/{df}/taxa/kraken2/{preproc}/{sample}_2.fq',
        o   = '{prefix}/{df}/taxa/kraken2/{preproc}/{sample}.krak',
        r   = '{prefix}/{df}/taxa/kraken2/{preproc}/{sample}.report',
    params: 
        out = '{prefix}/{df}/taxa/kraken2/{preproc}/{sample}#.fq'
    threads: 12
    conda: 'kraken2_env.yaml'
    shell: '''kraken2 --db /data5/bio/databases/kraken2/16S_Silva_20190418 --threads {threads} --fastq-input --gzip-compressed \
        --paired {input.r1} {input.r2} --classified-out {params.out} --use-mpa-style --report {output.r} > {output.o}'''