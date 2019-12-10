rule fastqc:
    input: wc_config['fastq_gz_file_wc']
    output: 
        zipped="{fs_prefix}/{df}/profile/{preproc}/{sample}/{sample}_{strand}_fastqc.zip"
    params: 
        out="{fs_prefix}/{df}/profile/{preproc}/{sample}/",
        zip_out="{fs_prefix}/{df}/reads/{preproc}/{sample}/profile/fastqc"
    log: "{fs_prefix}/{df}/profile/{preproc}/{sample}/{sample}_{strand}.log"
    threads: 6
    conda: 'env_v0.11.8.yaml'
    shell: ('''export PERL5LIB='';\nfastqc -t {threads} -o {params.zip_out} {input} >{log} 2>&1; \n
          unzip -o {output.zipped} -d {params.zip_out}''')
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