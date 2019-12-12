

rule kraken:
    input:
        r1 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz'
    output:
        report = 'datasets/{df}/taxa/{preproc}/kraken-{version}-ff/{db}/{sample}/report.tsv', # ff param is for FinfFungi
    params:
        classified = 'datasets/{df}/taxa/{preproc}/kraken-{version}-ff/{db}/{sample}/classified.tsv',
        unclassified = 'datasets/{df}/taxa/{preproc}/kraken-{version}-ff/{db}/{sample}/unclassified.tsv',
    threads:  12
    run:
        DB = config['kraken']['db'][str(wildcards.db)]
        KRAKEN = config["kraken"][str(wildcards.version)]
        shell ('''{KRAKEN} --db {DB} --threads {threads} \
         --preload \
         --output {output.report} \
         --classified-out {params.classified} \
         --fastq-input --gzip-compressed {input.r1} {input.r2}''')
