def get_samples_for_maxbin2_cov(wildcards):
    bam_wc = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/mapped.bb_stats'
    ss = wildcards.samples

    table_wc = '{prefix}/{df}/assembly/mh__{params}/{samples}/sample_set.tsv'
    r_wc_str = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz'
    
    table = pd.read_csv(table_wc.format(prefix = wildcards.prefix,
                                        df = wildcards.df,
                                        params = wildcards.params,
                                        samples = wildcards.samples),
                        sep = '\t')

    list_of_sample_profiles = []
    for s in table.to_dict(orient='records'):     
        list_of_sample_profiles.append(bam_wc.format(prefix=wildcards.prefix,
                                    bwa_params = wildcards.bwa_params,
                                    params = wildcards.params,
                                    mod =  wildcards.mod, 
                                    samples = wildcards.samples,
                                    df=s['df'], 
                                    preproc=s['preproc'],
                                    sample=s['fs_name']))

    print(list_of_sample_profiles)
    return list_of_sample_profiles 

rule prepare_cov_files:
    input: 
        fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.fa'),
        samples = get_samples_for_maxbin2_cov
    output:
        done = '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/covs/abund_list.txt'
    params:
        wd = '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/covs/'
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
        fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.fa'),
        abund_list = '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/covs/abund_list.txt'
    output:
        done       = '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/maxbin2.done'
    params:
        wd         = '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/bins/',
        prefix     = 'bin'
    log:             '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/log.txt'
    benchmark:       '{prefix}/{df}/maxbin2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/benchmark.txt'
    threads: 40
    conda: 'maxbin2_env.yaml'
    shell: ('''export PERL5LIB="/data6/bio/TFM/pipeline/.snakemake/conda/fcda6b9a/lib/site_perl/5.26.2";\n
        mkdir -p {params.wd}; \n
        cd {params.wd}; \n
        (run_MaxBin.pl -contig {input.fa} -out {params.prefix} -abund_list {input.abund_list} -thread {threads}) >{log} 2>&1; \n
            touch {output.done}''')