rule count:
    input: 
        r1 = "{fs_prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz",
    output: 
        r1 = "{fs_prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}.count"
    wildcard_constraints:    
        df="[\w\d_-]+"
    run:
        exec_loc = os.path.join(config['assnake_install_dir'], 'bin/scripts/count_bp.sh')
        shell('''{exec_loc} {input.r1} > {output.r1}''')

        #save_to_db(config['task_id'], rule, str(input), str(output.r1), 'RUN SUCCESSFUL')