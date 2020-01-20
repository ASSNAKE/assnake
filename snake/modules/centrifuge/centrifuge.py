CENTRIFUGE_INDEX = config["centrifuge"]["index"]
fna_db_dir = config['fna_db_dir']

rule centr_krona:
    input:classification = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_classification.tsv'
    output: krona = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krona.tsv'
    params: wd = 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/'
    # conda: 'env_v1.0.4_beta.yaml'
    run:
        shell('tail -n +2 {input.classification} | cut -f 1,3 > {output.krona}')
        shell('''cd {params.wd} \n /data5/bio/runs-fedorov/tools/krona/bin/ktImportTaxonomy ./{wildcards.sample}_krona.tsv -o krona.html''')


rule centrifuge:
    input: 
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        # params = "params/centr/{params}.json"
    output:
        report =         '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_report.tsv',
        classification = '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_classification.tsv',
        krak =           '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krak.tsv', 
    threads:  12
    log:       '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/log.txt'
    benchmark: '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/benchmark.txt'
    conda: 'env_v1.0.4_beta.yaml'
    shell: ('''centrifuge -k 1 --mm --min-hitlen 22 -f -x {CENTRIFUGE_INDEX} -1 {input.r1} -2 {input.r2} -p {threads} -S {output.classification} --report-file {output.report} -q >{log} 2>&1; \n
          centrifuge-kreport -x {CENTRIFUGE_INDEX} {output.classification} > {output.krak}''')

rule centrifuge_fasta:
    input: 
        fa = os.path.join(fna_db_dir,'{path}/{seq_set_id}.fa')
        # params = "params/centr/{params}.json"
    output:
        report =         os.path.join(fna_db_dir,'{path}/{seq_set_id}_centr__{params}/report.tsv'),
        classification = os.path.join(fna_db_dir,'{path}/{seq_set_id}_centr__{params}/classification.tsv'),
        krak =           os.path.join(fna_db_dir,'{path}/{seq_set_id}_centr__{params}/krak.tsv'), 
    threads:  12
    log:       os.path.join(fna_db_dir,'{path}/{seq_set_id}_centr__{params}/log.txt')
    benchmark: os.path.join(fna_db_dir,'{path}/{seq_set_id}_centr__{params}/benchmark.txt')
    conda: 'env_v1.0.4_beta.yaml'
    shell:
        ('''centrifuge -k 1 --mm --min-hitlen 120 -f -x {CENTRIFUGE_INDEX} -U {input.fa} -p {threads} -S {output.classification} --report-file {output.report} >{log} 2>&1; \n
         centrifuge-kreport -x {CENTRIFUGE_INDEX} {output.classification} > {output.krak}''')