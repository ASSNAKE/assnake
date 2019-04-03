CENTRIFUGE_FOLDER = config["centrifuge"]["bin"]
CENTRIFUGE_INDEX = config["centrifuge"]["index"]

rule centrifuge:
    input: 
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        # params = "params/centr/{params}.json"
    output:
        report =         '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_report.tsv',
        classification = '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_classification.tsv',
        krak =           '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krak.tsv', 
    threads:  6
    log:       '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/log.txt'
    benchmark: '{prefix}/{df}/taxa/{preproc}/centr__{params}/{sample}/benchmark.txt'
    run:
        shell('{CENTRIFUGE_FOLDER}centrifuge -k 1 --mm --min-hitlen 22 -f -x {CENTRIFUGE_INDEX} -1 {input.r1} -2 {input.r2} -p {threads} -S {output.classification} --report-file {output.report} -q >{log} 2>&1')
        shell('{CENTRIFUGE_FOLDER}centrifuge-kreport -x {CENTRIFUGE_INDEX} {output.classification} > {output.krak}')