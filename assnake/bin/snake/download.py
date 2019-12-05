rule fastq_dump:
    input: 
        sample     = '{prefix}/{df}/reads/sra/{sample}/{sample}.download'
    output:
        r1         = '{prefix}/{df}/reads/sra/{sample}/{sample}_R1.fastq.gz',
        r2         = '{prefix}/{df}/reads/sra/{sample}/{sample}_R2.fastq.gz'
    params: 
        out_folder = '{prefix}/{df}/reads/sra/{sample}/',
        r1         = '{prefix}/{df}/reads/sra/{sample}/{sample}_1.fastq.gz',
        r2         = '{prefix}/{df}/reads/sra/{sample}/{sample}_2.fastq.gz'
    log:             '{prefix}/{df}/reads/sra/{sample}/{sample}.download.log'
    run: 
        shell('/srv/common/bin/fastq-dump --split-files --gzip --dumpbase --skip-technical --outdir {params.out_folder} {wildcards.sample} >{log} 2>&1')
        shell('mv {params.r1} {output.r1}')
        shell('mv {params.r2} {output.r2}')