import shutil
import os
import pandas as pd
import assnake.api.dataset 

assembly_dir = config['assembly_dir']
assnake_db = config['assnake_db']
fna_db_dir = config['fna_db_dir']


def megahit_input_from_table(wildcards):
    """
    Reads table with samples, returns dict with reads
    """
    table_wc = '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/sample_set.tsv'
    r_wc_str = '{fs_prefix}/{df}/reads/{preproc}/{sample}_{strand}.fastq.gz'
    
    table = pd.read_csv(table_wc.format(fs_prefix = wildcards.fs_prefix,
                                        df = wildcards.df,
                                        params = wildcards.params,
                                        sample_set = wildcards.sample_set),
                        sep = '\t')
    print(table)
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


rule megahit_from_table:
    input:
        unpack(megahit_input_from_table),
        table      = '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/sample_set.tsv',
        params=os.path.join(config['assnake_db'], "params/megahit/{params}.json")
    output:
        out_fa     = '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/final_contigs.fa',
        params     = '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/params.json'
    params:
        out_folder = '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/assembly/'
    threads: 12
    log: '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/log.txt'
    conda: 'megahit_env_v1.2.9.yaml'
    wrapper: "file://"+os.path.join(config['assnake_install_dir'], 'modules/megahit/megahit_wrapper.py')

def get_ref(wildcards):
    fs_prefix = api.dataset.Dataset(wildcards.df).fs_prefix
    return '{fs_prefix}/{df}/assembly/mh__v1.2.9__{params}/{sample_set}/final_contigs.fa'.format(
            fs_prefix = fs_prefix, 
            df = wildcards.df,
            sample_set = wildcards.sample_set,
            params = wildcards.params)

rule refine_assemb_results_cross:
    input: 
        ref = '{fs_prefix}/{df}/assembly/{assembler}__{assembler_version}__{params}/{sample_set}/final_contigs.fa'
    output: 
        fa =         '{fs_prefix}/{df}/assembly/{assembler}__{assembler_version}__{params}/{sample_set}/final_contigs__{min_len}.fa',
        fai =        '{fs_prefix}/{df}/assembly/{assembler}__{assembler_version}__{params}/{sample_set}/final_contigs__{min_len}.fa.fai',
        dictionary = '{fs_prefix}/{df}/assembly/{assembler}__{assembler_version}__{params}/{sample_set}/final_contigs__{min_len}.fa.dict'
    log: 
        names = '{fs_prefix}/{df}/assembly/{assembler}__{assembler_version}__{params}/{sample_set}/final_contigs__{min_len}.txt',
        ll    = '{fs_prefix}/{df}/assembly/{assembler}__{assembler_version}__{params}/{sample_set}/final_contigs__{min_len}.log'
    wildcard_constraints:
        min_len="[\d_-]+"
    run:
        shell('echo -e "INFO: Filtering contigs < {wildcards.min_len}bp and simplifying names"')
        shell("{config[anvio.bin]}anvi-script-reformat-fasta {input.ref} -o {output.fa} --min-len {wildcards.min_len} --simplify-names --report {log.names} > {log.ll} 2>&1")
        shell('echo -e "INFO: Done filtering contigs < {wildcards.min_len} and simplifying names!"')
        shell('''echo -e "INFO: Creating fai and dict files for reference" \n
                {config[samtools.bin]} faidx {output.fa} \n
                {config[java.bin]} \
                    -jar {config[picard.jar]} CreateSequenceDictionary \
                        REFERENCE={output.fa} \
                        OUTPUT={output.dictionary} \n
                echo -e "\nINFO: Done creating fai and dict files for reference!"''')