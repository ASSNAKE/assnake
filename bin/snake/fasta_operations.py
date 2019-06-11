rule extract_sequence_names:
    input: os.path.join(fna_db_dir, '{path}/{seq_set_id}.fa')
    output: os.path.join(fna_db_dir, '{path}/{seq_set_id}.names')
    shell: ('''awk 'sub(/^>/, "")' {input} > {output}''')


rule filter_by_kumar:
    input: sam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.sam'
    output:
        bam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_kumar.bam'
    params:
        tmp = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_no_xa.tmp.bam',
        filt_sam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_filt.sam',
        no_xa = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_no_xa.sam',
    run: 
        shell('echo "{input.sam}"')
        shell('python {config[assnake_install_dir]}/assnake/bin/scripts/mgSNP_sam-filter.py -i {input.sam} -o {params.filt_sam} -m 97 -l 60')   
        shell('grep -v "XA:" {params.filt_sam} > {params.no_xa}')   
        shell('{config[samtools.bin]} view -bS {params.no_xa} -o {params.tmp}')
        shell('{config[samtools.bin]} sort {params.tmp} -o {output.bam}')
        shell('{config[samtools.bin]} index {output.bam}')
        shell('rm {params.tmp} {params.filt_sam} {params.no_xa}')