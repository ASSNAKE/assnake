FEATURE_COUNTS_BIN='/data6/bio/TFM/soft/bioinf/subread-1.6.3-source/bin/featureCounts'
nucl_dir= config['na_db_dir']
rule feature_counts:
    input: gff = nucl_dir+'{type}/{category}/{seq_object}/{seq_set_id}.gff',
           sam = 'datasets/{df}/mapped/{preproc}/bwa__{params}/{type}__{category}__{seq_object}/{seq_set_id}/mapped/{sample}.sam'
    output: fc = 'datasets/{df}/mapped/{preproc}/bwa__{params}/{type}__{category}__{seq_object}/{seq_set_id}/mapped/{sample}_feature_counts.txt'
    threads: 5
    run:
        shell('{FEATURE_COUNTS_BIN} -T {threads} -p -s 1 -t gene -g ID -a {input.gff} -o {output.fc} {input.sam}')