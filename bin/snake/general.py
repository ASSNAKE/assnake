fna_db_dir= config['fna_db_dir']



rule init_bam_mapped:
    input:
        sam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.sam'
    output:
        bam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.bam'
    params:
        tmp = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.tmp.bam'
    run: 
        shell('echo "{input.sam}"')   
        shell('{config[samtools.bin]} view -bS {input.sam} -o {params.tmp}')
        shell('{config[samtools.bin]} sort {params.tmp} -o {output.bam}')
        shell('{config[samtools.bin]} index {output.bam}')
        shell('rm {params.tmp}')

nucl_dir= config['fna_db_dir']
rule coverage:
    input:
        bam = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}/{id_seq_set}/{what}/{id_sample}.bam',
        ref = nucl_dir+"{type}/{fasta}.fa"
    output:
        cov = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}/{id_seq_set}/{what}/{id_sample}.cov'
    run:
        shell('/srv/common/bin/genomeCoverageBed -bga -ibam {input.bam} > {output.cov}')

def get_reference_fasta(wildcards):
    """
    Reconstructs reference fasta location
    """
    path = wildcards.path.replace('___', '/')
    
    fasta_wc_loc = os.path.join(fna_db_dir, '{path}/{seq_set_id}.fa')
    fasta = fasta_wc_loc.format(path=path, seq_set_id=wildcards.seq_set_id)

    return fasta

rule coverage_stats:
    input:
        bam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.bam',
        ref = get_reference_fasta
    output:
        stats = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.bb_stats'
    log: '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_bb_stats_log.txt'
    run:
        shell('''{config[java.bin]} -ea -Xmx35000m -cp /data6/bio/TFM/soft/bbmap/current/ jgi.CoveragePileup ref={input.ref} in={input.bam} out={output.stats} >{log} 2>&1 \n
        cat {log}''')
             
