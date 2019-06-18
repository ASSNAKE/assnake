rule fastqc:
    input: "{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz"
    output: 
        zipped="{prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}_fastqc.zip"
    params: 
        out="{prefix}/{df}/reads/{preproc}/{sample}/profile/",
        zip_out="{prefix}/{df}/reads/{preproc}/{sample}/profile/"
    log: "{prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}.log"
    threads: 6
    conda: 'env_v0.11.8.yaml'
    shell: ('''fastqc -t {threads} -o {params.out} {input} >{log} 2>&1; \n
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
    shell: ('''fastqc -t {threads} --nogroup -o {params.out} {input} >{log} 2>&1; \n
          unzip -o {output.zipped} -d {params.zip_out}''')