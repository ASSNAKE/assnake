rule just_a_test:
    input: "datasets/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz"
    output: "datasets/{df}/reads/{preproc}/{sample}/test/{sample}_{strand}_test.txt"
    run:
        shell('''echo "Some text here." > {output}''')
        save_to_db(config['task_id'], rule, str(input), str(output), 'RUN SUCCESSFUL')