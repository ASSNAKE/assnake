rule fastqc:
    input: wc_config['fastq_gz_file_wc']
    output: 
        zipped=wc_config['fastqc_zip_wc']
    params: 
        out="{fs_prefix}/{df}/profile/fastqc/{preproc}/{sample}/",
        zip_out="{fs_prefix}/{df}/profile/fastqc/{preproc}/{sample}/fastqc"
    log: "{fs_prefix}/{df}/profile/fastqc/{preproc}/{sample}/{sample}_{strand}.log"
    threads: 6
    conda: 'env_v0.11.8.yaml'
    shell: ('''export PERL5LIB='';\nfastqc -t {threads} -o {params.out} {input} >{log} 2>&1; \n
          unzip -o {output.zipped} -d {params.out}''')
        #save_to_db(config['task_id'], rule, str(input), str(output.zipped), 'RUN SUCCESSFUL')

rule fastqc_nogroup:
    input: "{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz"
    output: 
        zipped="{prefix}/{df}/reads/{preproc}/{sample}/profile/fastqc_nogroup/{sample}_{strand}_fastqc.zip"
    params: 
        out="{prefix}/{df}/reads/{preproc}/{sample}/profile/fastqc_nogroup/",
        zip_out="{prefix}/{df}/reads/{preproc}/{sample}/profile/fastqc_nogroup/"
    log: "{prefix}/{df}/reads/{preproc}/{sample}/profile/fastqc_nogroup/{sample}_{strand}.log"
    threads: 6
    conda: 'env_v0.11.8.yaml'
    shell: ('''export PERL5LIB='';\nfastqc -t {threads} --nogroup -o {params.out} {input} >{log} 2>&1; \n
          unzip -o {output.zipped} -d {params.zip_out}''')