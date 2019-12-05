def get_reference_fasta2(wildcards):
    """
    Reconstructs reference fasta location
    """
    path = wildcards.path.replace('___', '/')
    print('wtf')
    fasta_wc_loc = os.path.join(fna_db_dir, '{path}/{seq_set_id}.fa')
    fasta = fasta_wc_loc.format(path=wildcards.path, seq_set_id=wildcards.seq_set_id)
    print(fasta)
    # print(**wildcards)
    return fasta

rule index_bbmap:
    input:
        ref = get_reference_fasta2
    output:
        done = os.path.join(fna_db_dir, 'index/bbmap/{path}/{seq_set_id}/index.done')
    params:
        prefix    = os.path.join(fna_db_dir, 'index/bbmap/{path}/{seq_set_id}/')
    log:            os.path.join(fna_db_dir, 'index/bbmap/{path}/{seq_set_id}/log.txt')
    benchmark:      os.path.join(fna_db_dir, 'index/bbmap/{path}/{seq_set_id}/benchmark.txt')
    conda: 'bbmap_env.yaml'
    shell: ('''bbmap.sh ref={input.ref} path={params.prefix} -Xmx23g; \n
        touch {output.done}''')

rule coverage_stats:
    input:
        bam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.bam',
        ref = get_reference_fasta2
    output:
        stats = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.bb_stats'
    wildcard_constraints:    
        df="[\w\d_-]+",
        params="[\w\d_-]+"
    log: '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_bb_stats_log.txt'
    conda: 'bbmap_env.yaml'
    
    shell: ('''pileup.sh ref={input.ref} in={input.bam} out={output.stats} >{log} 2>&1 \n
    cat {log}''')