METAPHLAN2 = config['METAPHLAN2']
MPA_PKL = config['MPA_PKL']
BOWTIE2DB = config['BOWTIE2DB']
BOWTIE2 = config['bowtie2.bin']

CENTRIFUGE_FOLDER = config["centrifuge"]["bin"]
CENTRIFUGE_INDEX = config["centrifuge"]["index"]

rule metaphlan2:
    input: 
        r1 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        params = "params/mp2/{params}.json"
    output: 
        o = 'datasets/{df}/taxa/{preproc}/mp2__{params}/{sample}/{sample}.mp2'
    params:
        b =  'datasets/{df}/taxa/{preproc}/mp2__{params}/{sample}/{sample}.b2'
    threads: 6
    run:
        shell('''{METAPHLAN2} --mpa_pkl {MPA_PKL} --bowtie2db {BOWTIE2DB} --bowtie2_exe {BOWTIE2} {input.r1},{input.r2} --nproc {threads} --input_type fastq --bowtie2out {params.b}  > {output.o}''')
        if 'task_id' in config.keys():
            save_to_db(config['task_id'], 'mp2', str(input), str(output.o), 'RUN SUCCESSFUL')
               
rule centrifuge_contigs:
    input: 
        seq = 'data/ref/assemb/{assemb_name}.fa'
    output:
        report = 'data/taxa/assemb/{assemb_name}/centr_report.tsv',
        classification = 'data/taxa/assemb/{assemb_name}/centr_classification.tsv',
        krak = 'data/taxa/assemb/{assemb_name}/centr_krak.tsv',
        krona = 'data/taxa/assemb/{assemb_name}/centr_krona.tsv'
    params:
        wd = 'data/taxa/assemb/{assemb_name}/',
        krona_html = 'data/taxa/assemb/centr_krona'
    run:
        shell('{CENTRIFUGE_FOLDER}centrifuge -k 1 --min-hitlen 200 -f -x {CENTRIFUGE_INDEX} {input.seq} -S {output.classification} --report-file {output.report}')
        shell('{CENTRIFUGE_FOLDER}centrifuge-kreport -x {CENTRIFUGE_INDEX} {output.classification} > {output.krak}')
        shell('tail -n +2 {output.classification} | cut -f 1,3 > {output.krona}')
        shell('''cd {params.wd} \n
             /data5/bio/runs-fedorov/tools/krona/bin/ktImportTaxonomy ./centr_krona.tsv -o krona_centr.html''')
        
rule centrifuge:
    input: 
        r1 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        params = "params/centr/{params}.json"
    output:
        report =         'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_report.tsv',
        classification = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_classification.tsv',
        krak =           'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krak.tsv', 
    threads:  12
    log: 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/log.txt'
    benchmark: 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/time.txt'
    run:
        shell('{CENTRIFUGE_FOLDER}centrifuge -k 1 --mm --min-hitlen 22 -f -x {CENTRIFUGE_INDEX} -1 {input.r1} -2 {input.r2} -p {threads} -S {output.classification} --report-file {output.report} -q >{log} 2>&1')
        shell('{CENTRIFUGE_FOLDER}centrifuge-kreport -x {CENTRIFUGE_INDEX} {output.classification} > {output.krak}')
        
        
        
rule centr_krona:
    input:classification = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_classification.tsv'
    output: krona = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krona.tsv'
    params: wd = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/'
    run:
        shell('tail -n +2 {input.classification} | cut -f 1,3 > {output.krona}')
        shell('''cd {params.wd} \n /data5/bio/runs-fedorov/tools/krona/bin/ktImportTaxonomy ./{wildcards.sample}_krona.tsv -o krona.html''')




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