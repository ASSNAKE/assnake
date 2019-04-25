import os

nucl_dir= config['fna_db_dir']
fna_db_dir= config['fna_db_dir']

def get_reference_fasta(wildcards):
    """
    Reconstructs reference fasta location
    """
    path = wildcards.path.replace('___', '/')
    
    fasta_wc_loc = os.path.join(fna_db_dir, '{path}/{seq_set_id}.fa')
    fasta = fasta_wc_loc.format(path=path, seq_set_id=wildcards.seq_set_id)

    return fasta


rule create_seq_set_index_bwa:
    input:
        ref = get_reference_fasta
    output:
        ref_index = os.path.join(fna_db_dir, 'index/bwa/{path}/{seq_set_id}/index.sa')
    params:
        prefix    = os.path.join(fna_db_dir, 'index/bwa/{path}/{seq_set_id}/index')
    log:            os.path.join(fna_db_dir, 'index/bwa/{path}/{seq_set_id}/log.txt')
    benchmark:      os.path.join(fna_db_dir, 'index/bwa/{path}/{seq_set_id}/benchmark.txt')
    conda: 'env_0.7.17.yaml'
    shell: ('''echo -e "INFO: Creating BWA index for {input.ref}"; \n
         bwa index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1; \n
         echo -e "INFO: Finished creating BWA index for {input.ref}\n"''')
        

rule map_on_ref_bwa:
    input:
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        ref_index = os.path.join(fna_db_dir, 'index/bwa/{path}/{seq_set_id}/index.sa')
    output:
        sam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.sam'
    params:
        ind_prefix = os.path.join(fna_db_dir, 'index/bwa/{path}/{seq_set_id}/index')
    log:      '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/log.txt'
    benchmark:'{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/benchmark.txt'
    threads: 20
    conda: 'env_0.7.17.yaml'
    shell: ('(bwa mem -M -t {threads} {params.ind_prefix} {input.r1} {input.r2} | /srv/common/bin/samtools view -SF 4 -h > {output.sam}) >{log} 2>&1')

rule leave_aligned_bwa:
    input: 
        sam = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/{id_sample}.sam'
    output: 
        filtered = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/{id_sample}.aligned.sam'
    run:
        shell('samtools view -S -F 4 {input.sam} > {output.filtered}')
        

rule filter_mapping:
    input: 
        sam = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/mapped/{id_sample}.sam'
    output: 
        filtered = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/filtered/{id_sample}.sam'
    benchmark: 'time/filter-map/{df}/{preproc}/{type}/{id_seq_set}/{id_sample}.time'
    run:
        shell('{config[python.bin]} bin/scripts/mgSNP_sam-filter.py -i {input.sam} -o {output.filtered}.pre -l 25 -m 90')
        shell('grep -v "XA:" {output.filtered}.pre > {output.filtered}')
        shell('rm {output.filtered}.pre')
        

#maps reads from sample on provided contig and adds read groups.        
rule create_assemb_index_bwa:
    input:
        ref = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/final.contigs.fa'
    output:
        ind = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/index/bwa/index.sa'
    params:
        prefix = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/index/bwa/index'
    log: '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/index/bwa/log.txt'
    run:
        shell('echo -e "INFO: Creating BWA index for {input.ref}\n"')
        shell('{config[bwa.bin]} index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1')
        shell('echo -e "INFO: Finished creating BWA index for {input.ref}\n"')
        
rule map_on_contig_bwa:
    input:
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        ref_fasta = '{prefix}/{df}/assemb/mh__{params}/{preprocs}/{samples}/index/bwa/index.sa'
    output:
        sam = '{prefix}/{df}/mapped/{preproc}/bwa__{params}/assembRaw/{preprocs}__{samples}/mapped/{sample}.sam'
    log:      '{prefix}/{df}/mapped/{preproc}/bwa__{params}/assembRaw/{preprocs}__{samples}/mapped/{sample}.log'

    threads: 12
    run:
        ind_prefix = input.ref_fasta[0:-3]
        shell('({config[bwa.bin]} mem -M -t {threads} {ind_prefix} {input.r1} {input.r2} | /srv/common/bin/samtools view -SF 4 -h > {output.sam}) >{log} 2>&1')
        
rule compare_contigs:
    input: 
        ref_cont_ind = 'data/ref/index/bwa/assemb/{sample1}-{method}-{tool}/index.sa',
        cont1 = 'data/ref/assemb/{sample1}-{method}-{tool}.fa',
        cont2 = 'data/ref/assemb/{sample2}-{method}-{tool}.fa'
    output:
        comp_sam = 'data/ref/assemb/compare/{sample1}-{method}-{tool}_vs_{sample2}-{method}-{tool}.sam',
        comp_tsv = 'data/ref/assemb/compare/{sample1}-{method}-{tool}_vs_{sample2}-{method}-{tool}.tsv'
    params: 
        ref_pref = 'data/ref/index/bwa/assemb/{sample1}-{method}-{tool}/index',
        wd = 'data/ref/assemb/compare/'
    threads: 12
    log: 'log/bwa/map/compare/{sample1}-{method}-{tool}_vs_{sample2}-{method}-{tool}.log'
    run:
        shell('({config[bwa.bin]} mem -M -t {threads} {params.ref_pref} {input.cont2} > {output.comp_sam}) >{log} 2>&1')
        shell('{config[python.bin]} scripts/mgSNP_sam-filter.py -i {output.comp_sam} -o {params.wd}filtered.sam -m 95 -l 200')
        shell('grep -v "XA:" {params.wd}filtered.sam > {params.wd}filtered_no_xa.sam')
        shell('grep -v "@" {params.wd}filtered_no_xa.sam > {output.comp_tsv}')
        shell('sed -i "1i qname\tflag\trname\tpos\tmapq\tcigar\trnext\tpnext\tseq\ttlen\tqual\tNM\tMD\tAS\tXS\tSA\ttmp" {output.comp_tsv}')
        shell('rm {params.wd}filtered.sam {params.wd}filtered_no_xa.sam')

