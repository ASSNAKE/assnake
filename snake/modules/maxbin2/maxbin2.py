def get_samples_for_maxbin2_cov(wildcards):
    bam_wc = '{fs_prefix}/{df}/mapped/bwa__0.7.17__{bwa_params}/assembly/mh__v1.2.9__def/{ass_df}/{sample_set}/final_contigs__{mod}/{sample}/{preproc}/mapped.jgi_depth'
    ss = wildcards.sample_set

    table_wc = '{fs_prefix}/{df}/assembly/mh__v1.2.9__def/{sample_set}/sample_set.tsv'
    r_wc_str = wc_config['fastq_gz_file_wc']
    
    table = pd.read_csv(table_wc.format(fs_prefix = wildcards.fs_prefix,
                                        df = wildcards.ass_df,
                                        sample_set = wildcards.sample_set),
                        sep = '\t')

    list_of_sample_profiles = []
    for s in table.to_dict(orient='records'):     
        list_of_sample_profiles.append(bam_wc.format(fs_prefix=wildcards.fs_prefix,
                                    bwa_params = wildcards.bwa_params,
                                    # params = wildcards.params,
                                    mod =  wildcards.mod, 
                                    sample_set = wildcards.sample_set,
                                    ass_df = wildcards.ass_df,
                                    df=s['df'], 
                                    preproc=s['preproc'],
                                    sample=s['fs_name']))

    print(list_of_sample_profiles)
    return list_of_sample_profiles 

rule prepare_cov_files:
    input: 
        fa    = '{fs_prefix}/{ass_df}/assembly/mh__v1.2.9__def/{sample_set}/final_contigs__{mod}.fa',
        samples = get_samples_for_maxbin2_cov
    output:
        done = '{fs_prefix}/{df}/maxbin2/bwa__{version}__{bwa_params}/mh__v1.2.9__def/{ass_df}/{sample_set}/final_contigs__{mod}/covs/abund_list.txt'
    params:
        wd = '{fs_prefix}/{df}/maxbin2/bwa__{version}__{bwa_params}/mh__v1.2.9__def/{ass_df}/{sample_set}/final_contigs__{mod}/covs/'
    run:
        outs = [] 
        for s in input.samples:
            print(s.split('/')[-3])
            out = params.wd + s.split('/')[-3] + '.cov'
            outs.append(out)
            shell('''sed '1d' {s} | cut -f1,2 > {out}''')
        with open(output.done, 'w') as abund:
            abund.writelines('\n'.join(outs))
rule maxbin2:
    input:
        fa         = '{fs_prefix}/{df}/assembly/mh__v1.2.9__def/{sample_set}/final_contigs__{mod}.fa',
        abund_list = '{fs_prefix}/{df}/maxbin2/bwa__{version}__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/covs/abund_list.txt'
    output:
        done       = '{fs_prefix}/{df}/maxbin2/bwa__{version}__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/maxbin2.done'
    params:
        wd         = '{fs_prefix}/{df}/maxbin2/bwa__{version}__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/bins/',
        prefix     = 'bin'
    log:             '{fs_prefix}/{df}/maxbin2/bwa__{version}__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/log.txt'
    benchmark:       '{fs_prefix}/{df}/maxbin2/bwa__{version}__{params}/mh__v1.2.9__def/{df}/{sample_set}/final_contigs__{mod}/benchmark.txt'
    threads: 12
    conda: 'maxbin2_env.yaml'
    shell: ('''export PERL5LIB="/data6/bio/TFM/pipeline/.snake/conda/fcda6b9a/lib/site_perl/5.26.2";\n
        mkdir -p {params.wd}; \n
        cd {params.wd}; \n
        (run_MaxBin.pl -contig {input.fa} -out {params.prefix} -abund_list {input.abund_list} -thread {threads}) >{log} 2>&1; \n
            touch {output.done}''')