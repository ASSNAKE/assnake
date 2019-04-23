anvi_dir = 'data/anvi/'
ANVI = config['anvio.bin']
CENTRIFUGE_FOLDER = config["centrifuge"]["bin"]
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
    input: fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.fa')
    output: db_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.db')
    log: 
        # hmm = 'datasets/{df}/anvio/{type}/{ref}/db/hmm.log',
        gen =  os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr_db_gen_log.txt')
    threads: 10
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; anvi-gen-contigs-database -f {input.fa} -o {output.db_f} -n "{wildcards.samples}" >{log.gen} 2>&1''')
        # shell('''source /data6/bio/TFM/soft/miniconda3/bin/activate anvio5; anvi-run-hmms -c {output.db_f} -T {threads} > {log.hmm} 2>&1''')
#
# 
rule anvi_run_hmms:
    input: db_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.db')
    output: done = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr_hmms.done')
    log: 
        # hmm = 'datasets/{df}/anvio/{type}/{ref}/db/hmm.log',
        hmm =  os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr_db_hmm_log.txt')
    threads: 22
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; anvi-run-hmms -c {input.db_f} -T {threads} > {log.hmm} 2>&1''')     
        shell('touch {output.done}')    

rule anvi_cogs:
    input:db_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.db')
    output:cogs_done = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr_cogs.done')
    log: cogs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr_db_cogs_log.txt')
    threads: 22
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
    threads: 16
    run:
        shell("{ANVI}anvi-get-dna-sequences-for-gene-calls -c {input.db_f} -o {output.gc}")
        shell('{CENTRIFUGE_FOLDER}centrifuge -f -x {CENTRIFUGE_INDEX} {output.gc} -p {threads} -S {output.centr_hits} --report-file {output.centr_rep}')
        shell('{ANVI}anvi-import-taxonomy -c {input.db_f} -i {output.centr_rep} {output.centr_hits} -p centrifuge')
                
rule anvi_profile:
    input:
        bam = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly___mh__{params}___{dfs}___{samples}___{preprocs}/final_contigs__1000__no_hum_centr/{sample}/{preproc}/mapped.bam',
        contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.db')
        # centr_rep = anvi_dir+'db/{contig}-{method}-{tool}/centrifuge_report.tsv'
    output:
        aux    = '{prefix}/{df}/anvio/profile/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/{sample}___{preproc}/db/AUXILIARY-DATA.db',
        prof   = '{prefix}/{df}/anvio/profile/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/{sample}___{preproc}/db/PROFILE.db',
        log    = '{prefix}/{df}/anvio/profile/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/{sample}___{preproc}/db/RUNLOG.txt'
    params:
        dir_n =  '{prefix}/{df}/anvio/profile/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/{sample}___{preproc}/db/'
    log:         '{prefix}/{df}/anvio/profile/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/{sample}___{preproc}/log.txt'
    threads: 20
    run:
        shell('''
            source /data4/bio/fedorov/miniconda3/bin/activate anvio5; \n
            sample="{wildcards.sample}" \n
            sample=$(tr -d - <<< $sample) \n
                anvi-profile \
                    -i {input.bam} \
                    -c {input.contigs} \
                    --num-threads {threads} \
                    --sample-name $sample \
                    --output-dir {params.dir_n} \
                    --overwrite-output-destinations \
            >{log} 2>&1''')

                    #      --overwrite-output-destinations \


def get_samples_to_merge(wildcards):
    anvi_profile_wc = '{prefix}/{df}/anvio/profile/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/{sample}___{preproc}/db/PROFILE.db'
    
    samples = wildcards.samples.split(':')
    preproc = wildcards.preprocs

    list_of_sample_profiles = []
    for s in samples:
        list_of_sample_profiles.append(
            anvi_profile_wc.format(
                dfs = wildcards.dfs,
                samples = 'DFM_003_F1_S10:DFM_3F2_S63:D3T3_L1S1_B6_S106',
                preproc = wildcards.preprocs,
                preprocs = wildcards.preprocs,
                params = wildcards.params,
                bwa_params = wildcards.bwa_params,
                prefix = wildcards.prefix,
                df = wildcards.df,
                sample = s
            )
        )
    return list_of_sample_profiles
        
rule anvi_merge:
    input:
        contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.db'),
        samples = get_samples_to_merge
    output:        
        done = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/merge.done',
    log:
       '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/log.txt',
    params:
        out_dir = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/db',
    threads: 20
    run:
        shell('''source /data4/bio/fedorov/miniconda3/bin/activate anvio5; \n
            anvi-merge \
                {input.samples} \
                -o {params.out_dir} \
                -c {input.contigs} \
                --overwrite-output-destinations \
                --enforce-hierarchical-clustering \
                >{log} 2>&1''')
        shell('touch {output.done}')
        

        
rule anvi_get_hmm_seqs:
    input: db_f = anvi_dir+'db/{ref}/contigs.db'
    output: hmm_seqs = anvi_dir+'db/{ref}/hmm_seqs.fa'
    run: 
        shell('''{ANVI}anvi-get-sequences-for-hmm-hits -c {input.db_f} -o {output.hmm_seqs}''')
        
rule fuck:
    input:db_f = anvi_dir+'db/{ref}/contigs.db'
    output:func = anvi_dir+'db/{ref}/functions.txt'
    run: shell("{ANVI}anvi-export-functions -c {input.db_f} -o {output.func}")    
     