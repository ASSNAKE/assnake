def get_samples_for_maxbin2_cov(wildcards):
    bam_wc = '{prefix}/{df}/mapped/bwa__def1/assembly___mh__{params}___{dfs}___{samples}___{preprocs}/final_contigs__1000__no_hum_centr/{sample}/{preproc}/mapped.bb_stats'
    ss = wildcards.samples
    samples = wildcards.samples.split(':')
    preproc = wildcards.preprocs

    list_of_sample_profiles = []
    for s in samples:
        list_of_sample_profiles.append(
            bam_wc.format(
                prefix = wildcards.prefix,
                dfs = wildcards.dfs,
                samples = ss,
                preproc = wildcards.preprocs,
                preprocs = wildcards.preprocs,
                params = wildcards.params,
                bwa_params = wildcards.bwa_params,
                df = wildcards.df,
                sample = s
            )
        )
    return list_of_sample_profiles

def get_samples_for_maxbin2(wildcards):
    bam_wc = '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/covs/{sample}.cov'
    ss = wildcards.samples
    samples = wildcards.samples.split(':')
    preproc = wildcards.preprocs

    list_of_sample_profiles = []
    for s in samples:
        list_of_sample_profiles.append(
            bam_wc.format(
                prefix = wildcards.prefix,
                dfs = wildcards.dfs,
                samples = ss,
                preprocs = wildcards.preprocs,
                params = wildcards.params,
                bwa_params = wildcards.bwa_params,
                df = wildcards.df,
                sample = s
            )
        )
    return list_of_sample_profiles

rule prepare_cov_files:
    input: 
        fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.fa'),
        samples = get_samples_for_maxbin2_cov
    output:
        done = '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/covs/abund_list.txt'
    params:
        wd = '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/covs/'
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
        fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.fa'),
        abund_list = '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/covs/abund_list.txt'
    output:
        done       = '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/maxbin2.done'
    params:
        wd         = '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/bins/',
        prefix     = 'bin'
    log:             '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/log.txt'
    benchmark:       '{prefix}/{df}/maxbin2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/benchmark.txt'
    threads: 40
    conda: 'maxbin2_env.yaml'
    shell: ('''export PERL5LIB="/data6/bio/TFM/pipeline/.snakemake/conda/fcda6b9a/lib/site_perl/5.26.2";\n
        mkdir -p {params.wd}; \n
        cd {params.wd}; \n
        (run_MaxBin.pl -contig {input.fa} -out {params.prefix} -abund_list {input.abund_list} -thread {threads}) >{log} 2>&1; \n
            touch {output.done}''')