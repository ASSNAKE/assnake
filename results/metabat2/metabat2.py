def get_samples_for_metabat(wildcards):
    bam_wc = '{prefix}/{df}/mapped/bwa__def1/assembly___mh__{params}___{dfs}___{samples}___{preprocs}/final_contigs__1000__no_hum_centr/{sample}/{preproc}/mapped.bam'
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
    print(list_of_sample_profiles)
    return list_of_sample_profiles

rule metabat2:
    input:
        fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.fa'),
        samples = get_samples_for_metabat
    output:
        done = '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/done'
    params:
        wd = '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/'
    log:      '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/log.txt'
    benchmark:'{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/benchmark.txt'
    threads: 20
    conda: 'metabat2_env.yaml'
    shell: ('''
            cd {params.wd}; \n
            (runMetaBat.sh {input.fa} {input.samples}) >{log} 2>&1; \n
            touch {output.done}''')