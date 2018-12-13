anvi_dir = 'data/anvi/'
ANVI = config['anvio.bin']
CENTRIFUGE_FOLDER = config["centrifuge"]["bin"]
CENTRIFUGE_INDEX = config["centrifuge"]["index"]

rule anvi_gen_cont_db:
    input: fa = 'data/ref/{type}/{ref}.fa'
    output: db_f = 'datasets/{df}/anvio/{type}/{ref}/db/contigs.db'
    log: 
        hmm = 'datasets/{df}/anvio/{type}/{ref}/db/hmm.log',
        gen = 'datasets/{df}/anvio/{type}/{ref}/db/gen.log'
    threads: 10
    run:
        shell('''source /data6/bio/TFM/soft/miniconda3/bin/activate anvio5; anvi-gen-contigs-database -f {input.fa} -o {output.db_f} -n "{wildcards.ref}" >{log.gen} 2>&1''')
        shell('''source /data6/bio/TFM/soft/miniconda3/bin/activate anvio5; anvi-run-hmms -c {output.db_f} -T {threads} > {log.hmm} 2>&1''')
#         

rule anvi_cogs:
    input:db_f = 'datasets/{df}/anvio/{type}/{ref}/db/contigs.db'
    output:cogs_done = 'datasets/{df}/anvio/{type}/{ref}/db/cogs.done'
    log: cogs = 'datasets/{df}/anvio/{type}/{ref}/db/cogs.log'
    threads: 16
    run:
        shell('''{ANVI}anvi-run-ncbi-cogs --cog-data-dir {config[cogs.dir]} -T {threads} -c {input.db_f} >{log.cogs} 2>&1''')
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
        bam = 'data/mapped/bwa/assemb/{sample}_vs_{contig}-{method}-{tool}.bam',
        contigs = anvi_dir+'db/{contig}-{method}-{tool}/contigs.db',
        centr_rep = anvi_dir+'db/{contig}-{method}-{tool}/centrifuge_report.tsv'
    output:
        aux=anvi_dir+'profile/{contig}-{method}-{tool}/{sample}/AUXILIARY-DATA.db',
        prof=anvi_dir+'profile/{contig}-{method}-{tool}/{sample}/PROFILE.db',
        log=anvi_dir+'profile/{contig}-{method}-{tool}/{sample}/RUNLOG.txt'
    params:
        dir_n = anvi_dir+'profile/{contig}-{method}-{tool}/{sample}/'
    log: 'log/anvi/profile/{contig}-{method}-{tool}-{sample}.log'
    threads: 5
    run:
        shell('''
            sample="{wildcards.sample}" \n
            sample=$(tr -d - <<< $sample) \n
                {ANVI}anvi-profile \
                    -i {input.bam} \
                    -c {input.contigs} \
                    --num-threads {threads} \
                    --sample-name $sample \
                    --output-dir {params.dir_n} \
                    --overwrite-output-destinations \
            >{log} 2>&1''')

def get_samples_to_merge(wildcards):
    anvi_dir = 'data/anvi/'
    prefix = anvi_dir+'profile/'+wildcards.contig+'-'+wildcards.method+'-'+wildcards.tool+'/'
    samples = wildcards.samples.split(':')
    list_of_sample_profiles = []
    for s in samples:
        list_of_sample_profiles.append(prefix + s + '/PROFILE.db')
    return list_of_sample_profiles
        
rule anvi_merge:
    input:
        contigs = anvi_dir+'db/{contig}-{method}-{tool}/contigs.db',
        samples = get_samples_to_merge
    output:        
        done = anvi_dir+'profile/{contig}-{method}-{tool}/SAMPLES_MERGED/{samples}.done'
    log:
        'log/anvi/profile/{contig}-{method}-{tool}-{samples}-MERGED.log'
    params:
        out_dir = anvi_dir + 'profile/{contig}-{method}-{tool}/SAMPLES_MERGED/'
    run:
        shell('''
            {ANVI}anvi-merge \
                {input.samples} \
                -o {params.out_dir} \
                -c {input.contigs} \
                --overwrite-output-destinations \
            >{log} 2>&1
        ''')
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
     