rule bbmap_minlen:
    input:
        r1 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}/{sample}_R2.fastq.gz',
    output:
        r1 = '{fs_prefix}/{df}/reads/{preproc}__minlen{minlen}/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}__minlen{minlen}/{sample}_R2.fastq.gz',
    conda: '../bbmap/bbmap_env.yaml'
    shell: ('''reformat.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} minlength={wildcards.minlen}''')
