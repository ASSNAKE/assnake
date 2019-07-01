rule mothur_make_contigs:
    input: 
        samples = '{prefix}/{df}/mothur/{mothur_set}/stability.files'
    output: done = '{prefix}/{df}/mothur/{mothur_set}/merging.done'
    params: wd = '{prefix}/{df}/mothur/{mothur_set}/'
    log: '{prefix}/{df}/mothur/{mothur_set}/merging.log'
    threads: 20
    conda: 'mothur_env.yaml'
    shell: ('''cd {params.wd};\n
        mothur "#make.contigs(file={input.samples}, processors={threads})"  >{log} 2>&1;\n
        touch {output.done}''') 
