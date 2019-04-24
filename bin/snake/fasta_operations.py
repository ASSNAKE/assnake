rule extract_sequence_names:
    input: os.path.join(fna_db_dir, '{path}/{seq_set_id}.fa')
    output: os.path.join(fna_db_dir, '{path}/{seq_set_id}.names')
    shell: ('''awk 'sub(/^>/, "")' {input} > {output}''')