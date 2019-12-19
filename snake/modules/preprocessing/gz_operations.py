rule ungz_file:
    input: wc_config['fastq_gz_file_wc'] 
    output: wc_config['fastq_file']
    shell: ('''gunzip -c {input} > {output}''')

rule md5_sum_file:
    input: wc_config['fastq_gz_file_wc'] 
    output: wc_config['md5_fastq_gz']
    shell: ('''md5sum {input} > {output}''')