rule index_bbmap:
    input:
        ref = get_reference_fasta
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
        ref = get_reference_fasta
    output:
        stats = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.bb_stats'
    log: '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_bb_stats_log.txt'
    conda: 'bbmap_env.yaml'
    
    shell: ('''pileup.sh ref={input.ref} in={input.bam} out={output.stats} >{log} 2>&1 \n
    cat {log}''')