rule trim_galore:
    input: 
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
    output:
        r1 = '{prefix}/{df}/reads/{preproc}__trimgalore/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}__trimgalore/{sample}/{sample}_R2.fastq.gz',
        # dumb = '{prefix}/{df}/reads/{preproc}__trimgalore/{sample}/{sample}.done'
    params:
        wd = '{prefix}/{df}/reads/{preproc}__trimgalore/{sample}/',
        r1 = '{prefix}/{df}/reads/{preproc}__trimgalore/{sample}/{sample}_R1_val_1.fq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}__trimgalore/{sample}/{sample}_R2_val_2.fq.gz'
    threads: 4
    conda: 'env_0.6.2.yaml'
    shell: ('''trim_galore --nextera -o {params.wd} -j {threads} --paired {input.r1} {input.r2};\n 
        mv {params.r1} {output.r1}; mv {params.r2} {output.r2}''')