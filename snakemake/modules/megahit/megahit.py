import shutil
import os

def megahit_input(wildcards):
    """
    Generates input files for megahit
    {prefix}/{df}/assemb/mh__{params}/raw+imp/t:d+f:w:q/final.contigs.fa
    
    :param wildcards:
    :return:
    """
    print(wildcards.prefix)
    print(wildcards.df)

    r_wc_str = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz'
    
    # Reconstruct dict
    preprocs = wildcards.preprocs.split('+')
    samples_by_preproc = wildcards.samples.split('+')
    
    rr1 = []
    rr2 = []
    for i, preproc in enumerate(preprocs):
        samples = samples_by_preproc[i].split(':')
        for s in samples:
            rr1.append(r_wc_str.format(prefix=wildcards.prefix, 
                                       df = wildcards.df, 
                                       preproc=preproc,
                                       sample=s,
                                       strand='R1'))
            rr2.append(r_wc_str.format(prefix=wildcards.prefix, 
                                       df = wildcards.df, 
                                       preproc=preproc,
                                       sample=s,
                                       strand='R2'))
        
    print(rr1)
    return {'F': rr1, 'R': rr2}
        

rule run_megahit:
    input: 
        unpack(megahit_input)
    output: 
        out_fa = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/final_contigs.fa'
    params:
        out_folder = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/assemb/'
    threads: 12
    log: '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/log.txt'
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
        ref = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/final_contigs.fa'
    output: 
        fa =         '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/final_contigs_refined.fa',
        fai =        '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/final_contigs_refined.fa.fai',
        dictionary = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/final_contigs_refined.dict'
    log: 
        names = "{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/name_conversions.txt",
        ll = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/reformat_fasta.log'
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