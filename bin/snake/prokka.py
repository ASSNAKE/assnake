rule prokka:
    input: 'data/ref/{type}/{ref}/{ref}.fa'
    output: 'data/ref/{type}/{ref}/prokka'
    threads: 24
    shell: """ source /data6/bio/TFM/soft/miniconda3/bin/activate prokka; --cpus {threads} --outdir {output} --force --prefix {wildcards.ref} --locustag {wildcards.ref} {output} """