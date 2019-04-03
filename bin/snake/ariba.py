rule ariba_run:
    input:
        r1 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        db = '/data5/bio/databases/ariba/card/out.card.prepareref'
    output:
        res = 'datasets/{df}/ariba/card/{sample}/{preproc}/res.done',
        ar = 'datasets/{df}/ariba/card/{sample}/{preproc}/ariba.txt',
    conda: "/data6/bio/TFM/pipeline/envs/ariba.yml"
    shell: ('''ariba run {input.db} {input.r1} {input.r2} {output.ar} \n
     touch {output.res}''')
