import assnake.api.loaders

anvi_dir = 'data/anvi/'
# ANVI = config['anvio.bin']
# CENTRIFUGE_FOLDER = config["centrifuge"]["bin"]
CENTRIFUGE_INDEX = config["centrifuge"]["index"]
fna_db_dir= config['fna_db_dir']

clean_script = os.path.join( config['assnake_install_dir'], 'bin/scripts/filter_contigs_using_centrifuge.py')
human_contigs_list_script = os.path.join( config['assnake_install_dir'], 'bin/scripts/human_contigs_list.py')

rule clean_contigs_from_human:
    input:
        ll = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__human_contigs.list'),
        fa = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000.fa'),
    output: fa = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.fa')
    conda: 'env.yaml'
    shell: ('''seqkit grep -v -n -f {input.ll} {input.fa} > {output.fa}''')

rule centr_human_contigs:
    input:
        fa = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000.fa'),
        centr = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__centr__def_classification.tsv'),
    output: ll = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__human_contigs.list')
    shell: ('''python3 {human_contigs_list_script} --classification {input.centr} --out {output.ll}''')

rule anvi_gen_cont_db:
    input: fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.fa')
    output: db_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db')
    log: 
        # hmm = 'datasets/{df}/anvio/{type}/{ref}/db/hmm.log',
        gen =  os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db_gen.log')
    threads: 24
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; anvi-gen-contigs-database -f {input.fa} -o {output.db_f} -n "{wildcards.samples}" >{log.gen} 2>&1''')
        # shell('''source /data6/bio/TFM/soft/miniconda3/bin/activate anvio5; anvi-run-hmms -c {output.db_f} -T {threads} > {log.hmm} 2>&1''')
#
# 
rule anvi_run_hmms:
    input: db_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db')
    output: done = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}_hmms.done')
    log: 
        # hmm = 'datasets/{df}/anvio/{type}/{ref}/db/hmm.log',
        hmm =  os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}_hmm.log')
    threads: 24
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; anvi-run-hmms -c {input.db_f} -T {threads} > {log.hmm} 2>&1''')     
        shell('touch {output.done}')    

rule anvi_cogs:
    input:
        db_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db'),
        hmms = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}_hmms.done')
    output:cogs_done = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}_cogs.done')
    log: cogs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}_cogs.log')
    threads: 24
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5;  anvi-run-ncbi-cogs --cog-data-dir {config[cogs_dir]} -T {threads} -c {input.db_f} >{log.cogs} 2>&1''')

        shell('touch {output.cogs_done}')
        
rule anvi_centrifuge:
    input:
        db_f = anvi_dir+'db/{ref}/contigs.db'
    output:
        gc = anvi_dir+'db/{ref}/gene-calls.fa',
        centr_hits = anvi_dir+'db/{ref}/centrifuge_hits.tsv',
        centr_rep = anvi_dir+'db/{ref}/centrifuge_report.tsv'
    threads: 24
    run:
        shell("{ANVI}anvi-get-dna-sequences-for-gene-calls -c {input.db_f} -o {output.gc}")
        shell('{CENTRIFUGE_FOLDER}centrifuge -f -x {CENTRIFUGE_INDEX} {output.gc} -p {threads} -S {output.centr_hits} --report-file {output.centr_rep}')
        shell('{ANVI}anvi-import-taxonomy -c {input.db_f} -i {output.centr_rep} {output.centr_hits} -p centrifuge')
                
rule anvi_profile:
    input:
        bam = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/mapped.bam',
        contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db')
        # centr_rep = anvi_dir+'db/{contig}-{method}-{tool}/centrifuge_report.tsv'
    output:
        aux    = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/db/AUXILIARY-DATA.db',
        prof   = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/db/PROFILE.db',
        log    = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/db/RUNLOG.txt'
    params:
        dir_n =  '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/db/'
    log:         '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/db/log.txt'
    threads: 20
    run:
        sample_name = wildcards.sample
        if str.isdigit(sample_name[0]):
            sample_name = 's' + sample_name 
        shell('''
            source /data4/bio/fedorov/miniconda3/bin/activate anvio5; \n
            sample="{sample_name}"; \n
            sample=$(tr -d - <<< $sample); \n
            anvi-profile \
                -i {input.bam} \
                -c {input.contigs} \
                --num-threads {threads} \
                --sample-name $sample \
                --output-dir {params.dir_n} \
                --overwrite-output-destinations''')
#>{log} 2>&1
                    #      --overwrite-output-destinations \


def get_samples_to_merge(wildcards):
    anvi_profile_wc = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/db/PROFILE.db'
    ss = wildcards.samples

    table_wc = '{prefix}/{df}/assembly/mh__{params}/{samples}/sample_set.tsv'
    
    table = pd.read_csv(table_wc.format(prefix = wildcards.prefix,
                                        df = wildcards.df,
                                        params = wildcards.params,
                                        samples = wildcards.samples),
                        sep = '\t')

    list_of_sample_profiles = []
    for s in table.to_dict(orient='records'):     
        list_of_sample_profiles.append(anvi_profile_wc.format(prefix=wildcards.prefix,
                                    bwa_params = wildcards.bwa_params,
                                    params = wildcards.params,
                                    mod =  wildcards.mod, 
                                    samples = wildcards.samples,
                                    df=s['df'], 
                                    preproc=s['preproc'],
                                    sample=s['fs_name']))

    print(list_of_sample_profiles)
    return list_of_sample_profiles
        
rule anvi_merge:
    input:
        contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db'),
        # cogs_done = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}_cogs.done'),
        samples = get_samples_to_merge
    output:        
        done   = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}/mh__{params}/{df}/{samples}/{mod}/MERGED/merge.done',
        merged = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}/mh__{params}/{df}/{samples}/{mod}/MERGED/db/PROFILE.db',
    log: '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}/mh__{params}/{df}/{samples}/{mod}/MERGED/log.txt',
    params:
        out_dir = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}/mh__{params}/{df}/{samples}/{mod}/MERGED/db',
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; \n
            anvi-merge \
                {input.samples} \
                -o {params.out_dir} \
                -c {input.contigs} \
                --overwrite-output-destinations \
                --skip-concoct-binning \
                >{log} 2>&1''')
        shell('touch {output.done}')
        

# anvi-import-collection -p ./db/PROFILE.db -c /data5/bio/databases/fna/assembly/mh__def/FHM/D2/final_contigs__1000.db -C METABAT2 --contigs-mode /data10/bio/metagenomics/FHM/metabat2/bwa__def/mh__def/D2/1000/external_binning_of_contigs.tsv

def get_bins(wildcards):
    bin_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/conocot_anvio5_def/bin_by_bin')

def get_merged_profile(wildcards):
    dfs = api.loaders.load_dfs_from_db('')
    df_info = dfs[wildcards.df]

    merged_wc = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}/mh__{params}/{df}/{samples}/{mod}/MERGED/db/PROFILE.db'
    merged = merged_wc.format(
        prefix = df_info['fs_prefix'],
        df = wildcards.df,
        bwa_params = 'def',
        params = wildcards.params,
        samples = wildcards.samples,
        mod = wildcards.mod
    )
    return merged

rule anvi_summarize:
    input:
        contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.db'),
        merged = get_merged_profile
    output:        
        bins_summary = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/bins_summary.txt'),
    params:
        wd = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/'),
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; \n
                rm -rf {params.wd}; \n
                anvi-summarize \
                 -p {input.merged} \
                 -c {input.contigs} \
                 -C {wildcards.collection} \
                 -o {params.wd}''')
        
rule export_cov:
    input:
        contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.db'),
        merged = '/data6/bio/TFM/pipeline/datasets/{dfs}/anvio/merged_profiles/bwa__def1___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/db/PROFILE.db',
    output: '/data6/bio/TFM/pipeline/datasets/{dfs}/anvio/merged_profiles/bwa__def1___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/cov/contigs-COVs.txt'
    params: '/data6/bio/TFM/pipeline/datasets/{dfs}/anvio/merged_profiles/bwa__def1___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/cov',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; \n
        anvi-export-splits-and-coverages -p {input.merged} -c {input.contigs} -o {params} -O contigs --report-contigs''')
