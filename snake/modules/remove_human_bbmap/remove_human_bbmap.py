rule remove_human_bbmap:
    input:
        r1 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R2.fastq.gz',
        index = os.path.join(fna_db_dir, 'index/bbmap/human/hg19_bbmask/index.done')
    output:
        clean = '{fs_prefix}/{df}/reads/{preproc}__rmhum_bbmap/{sample}.fastq.gz',
    params:
        prefix    = os.path.join(fna_db_dir, 'index/bbmap/human/hg19_bbmask/')
    threads: 8
    conda: '../bbmap/bbmap_env.yaml'
    shell: ('''bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path={params.prefix} qtrim=rl trimq=10 untrim -Xmx23g in={input.r1} in2={input.r2} outu={output.clean} -t={threads}''')

#reformat.sh in=clean.fq out1=clean_R1.fq out2=clean_R2.fq

rule reformat_bbmap_hum:
    input:
        clean = '{fs_prefix}/{df}/reads/{preproc}__rmhum_bbmap/{sample}.fastq.gz',
    output:
        r1 = '{fs_prefix}/{df}/reads/{preproc}__rmhum_bbmap/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}__rmhum_bbmap/{sample}_R2.fastq.gz',
    conda: '../bbmap/bbmap_env.yaml'
    shell: ('''reformat.sh in={input.clean} out1={output.r1} out2={output.r2}''')
