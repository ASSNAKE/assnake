import shutil
import os

def megahit_input(wildcards):
    r_wc_str = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz'
    # Reconstruct dict
    dfs = wildcards.dfs.split('+')
    preprocs_by_df = wildcards.preprocs.split('+')
    samples_by_df = wildcards.samples.split('+')
    rrr = []
    
    rr1 = []
    rr2 = []
    for i, df in enumerate(dfs):
        rrr.append({df:[]})
        pr = preprocs_by_df[i].split('=')
        samps_by_df_preproc = samples_by_df[i].split('=')
        for j, p in enumerate(pr):
            rrr[i][df].append({p:[]})
            samps = samps_by_df_preproc[j].split(':')
            for s in samps:
                rrr[i][df][j][p].append(s)
                rr1.append(r_wc_str.format(df=df, preproc=p,sample=s,strand='R1'))
                rr2.append(r_wc_str.format(df=df, preproc=p,sample=s,strand='R2'))
    return {'F': rr1, 'R': rr2}
        

rule run_megahit:
    input: 
        unpack(megahit_input)
    output: 
        out_fa = 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/final.contigs.fa'
    params:
        out_folder = 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/assemb/'
    threads: 12
    log: 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/log.txt'
    run:
        reads1 = ",".join(input.F)
        reads2 = ",".join(input.R)
        print(reads1)
        if os.path.exists(params.out_folder) and os.path.isdir(params.out_folder):
            shutil.rmtree(params.out_folder)
            #os.makedirs(params.out_folder)
        shell('{config[megahit.bin]} -1 {reads1} -2 {reads2} --min-contig-len 850 -o {params.out_folder} -t {threads}  >{log} 2>&1')
        fc_loc = params.out_folder+'final.contigs.fa'
        shell('cp {fc_loc} {output.out_fa}')

rule refine_assemb_results:
    input: 
        ref = 'datasets/{df}/assemb/{preproc}/mh/{sample}/final.contigs.fa'
    output: 
        fa =         'data/ref/assemb-{df}-{preproc}-{tool}/{sample}.fa',
        fai =        'data/ref/assemb-{df}-{preproc}-{tool}/{sample}.fa.fai',
        dictionary = 'data/ref/assemb-{df}-{preproc}-{tool}/{sample}.dict'
    log: 
        names = "datasets/{df}/assemb/{preproc}/{tool}/{sample}/name_conversions.txt",
        ll = 'datasets/{df}/assemb/{preproc}/{tool}/{sample}/reformat_fasta.log'
    run:
        shell('echo -e "INFO: Filtering contigs < {config[min_contig_size]}bp and simplifying names"')
        shell("{config[anvi.bin]}anvi-script-reformat-fasta {input.ref} -o {output.fa} --min-len {config[min_contig_size]} --simplify-names --report {log.names} > {log.ll} 2>&1")
        shell('echo -e "INFO: Done filtering contigs < {config[min_contig_size]} and simplifying names!"')
        shell('''echo -e "INFO: Creating fai and dict files for reference" \n
                {config[samtools.bin]} faidx {output.fa} \n
                {config[java.bin]} \
                    -jar {config[picard.jar]} CreateSequenceDictionary \
                        REFERENCE={output.fa} \
                        OUTPUT={output.dictionary} \n
                echo -e "\nINFO: Done creating fai and dict files for reference!"''')