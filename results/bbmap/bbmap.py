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