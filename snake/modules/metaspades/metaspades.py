
def metaspades_input_from_table(wildcards):
    """
    Reads table with samples, returns dict with reads
    """
    table_wc = '{fs_prefix}/{df}/assembly/{sample_set}/sample_set.tsv'
    r_wc_str = '{fs_prefix}/{df}/reads/{preproc}/{sample}_{strand}.fastq.gz'
    
    table = pd.read_csv(table_wc.format(fs_prefix = wildcards.fs_prefix,
                                        df = wildcards.df,
                                        sample_set = wildcards.sample_set), sep = '\t')
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

rule prepare_merged_library:
    input: 
        unpack(metaspades_input_from_table),
        sample_set='{fs_prefix}/{df}/assembly/{sample_set}/sample_set.tsv',
    output: 
        r1 = '{fs_prefix}/{df}/assembly/{sample_set}/lib_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/assembly/{sample_set}/lib_R2.fastq.gz'
    wildcard_constraints:    
        df="[\w\d_-]+",
        preproc="[\w\d_-]+",
    shell: ('echo "{input.F}"; echo "{input.R}"; zcat {input.F} | gzip -c > {output.r1}; zcat {input.R} | gzip -c > {output.r2}')

rule metaspades:
    input:
        r1 = '{fs_prefix}/{df}/assembly/{sample_set}/lib_R1.fastq.gz',
        r2 = '{fs_prefix}/{df}/assembly/{sample_set}/lib_R2.fastq.gz',
        table      = '{fs_prefix}/{df}/assembly/{sample_set}/sample_set.tsv',
        params=os.path.join(config['assnake_db'], "params/metaspades/v3.14.0/{params}.json")
    output: 
        done = '{fs_prefix}/{df}/assembly/{sample_set}/metaspades__v3.14.0__{params}/done.txt'
    params: 
        out_folder = '{fs_prefix}/{df}/assembly/{sample_set}/metaspades__v3.14.0__{params}/assembly' 
    log: '{fs_prefix}/{df}/assembly/{sample_set}/metaspades__v3.14.0__{params}/log.txt'
    conda: 'metaspades_env.yaml'
    threads: 20
    wrapper: "file://"+os.path.join(config['assnake_install_dir'], 'modules/metaspades/metaspades_wrapper.py')