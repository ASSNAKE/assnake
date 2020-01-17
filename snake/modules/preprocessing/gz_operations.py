rule ungz_file:
    input: wc_config['fastq_gz_file_wc'] 
    output: wc_config['fastq_file']
    shell: ('''gunzip -c {input} > {output}''')

rule md5_sum_file:
    input: wc_config['fastq_gz_file_wc'] 
    output: wc_config['md5_fastq_gz']
    shell: ('''md5sum {input} > {output}''')

def get_files(wildcards):

    table_wc = '{fs_prefix}/{df}/reads/{preproc}__merge/{sample}_merge.tsv'
    r_wc_str = '{fs_prefix}/{df}/reads/{preproc}/{sample}_{strand}.fastq.gz'

    table = pd.read_csv(table_wc.format(fs_prefix = wildcards.fs_prefix,
                                        df = wildcards.df,
                                        preproc = wildcards.preproc,
                                        sample = wildcards.sample),
                        sep = '\t')
    rr1 = []
    rr2 = []                        
    for s in table.to_dict(orient='records'):     
        rr1.append(r_wc_str.format(fs_prefix=wildcards.fs_prefix,
                                    df=s['df'], 
                                    preproc=s['preproc'],
                                    sample=s['fs_name'],
                                    strand='R1'))
        rr2.append(r_wc_str.format(fs_prefix=wildcards.fs_prefix,
                                    df=s['df'], 
                                    preproc=s['preproc'],
                                    sample=s['fs_name'],
                                    strand='R2'))
    return {'F': rr1, 'R': rr2}

rule merge_fastq_gz_files:
    input: 
        unpack(get_files),
        sample_set='{fs_prefix}/{df}/reads/{preproc}__merge/{sample}_merge.tsv',
    output: 
        r1 = '{fs_prefix}/{df}/reads/{preproc}__merge/{sample}_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/reads/{preproc}__merge/{sample}_R2.fastq.gz'
    wildcard_constraints:    
        df="[\w\d_-]+",
        preproc="[\w\d_-]+",
    shell: ('echo "{input.F}"; echo "{input.R}"; zcat {input.F} | gzip -c > {output.r1}; zcat {input.R} | gzip -c > {output.r2}')