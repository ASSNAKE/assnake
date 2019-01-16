nucl_dir= config['na_db_dir']

# {type}/{category}/{seq_object}/{seq_set_id}.fa
# type = handplaced | imported | ??
# category = anything
rule create_seq_set_index_bwa:
    input:
        ref = nucl_dir+'{type}/{category}/{seq_object}/{seq_set_id}.fa'
    output:
        ref_index = nucl_dir+'index/bwa/{type}/{category}/{seq_object}/{seq_set_id}/index.sa'
    params:
        prefix = nucl_dir+'index/bwa/{type}/{category}/{seq_object}/{seq_set_id}/index'
    log: nucl_dir+'index/bwa/{type}/{category}/{seq_object}/{seq_set_id}/log.txt'
    benchmark: nucl_dir+'index/bwa/{type}/{category}/{seq_object}/{seq_set_id}/time.txt'
    run:
        shell('echo -e "INFO: Creating BWA index for {input.ref}\n"')
        shell('{config[bwa.bin]} index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1')
        shell('echo -e "INFO: Finished creating BWA index for {input.ref}\n"')
        
rule create_assemb_index_bwa:
    input:
        ref = 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/final.contigs.fa'
    output:
        ind = 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/index/bwa/index.sa'
    params:
        prefix = 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/index/bwa/index'
    log: 'processing/assemb/mh__{params}/{dfs}/{preprocs}/{samples}/index/bwa/log.txt'
    
    #benchmark: 'time/bwa/index/{type}/{contig}.time.txt'
    run:
        shell('echo -e "INFO: Creating BWA index for {input.ref}\n"')
        shell('{config[bwa.bin]} index -p {params.prefix} -a bwtsw {input.ref} > {log} 2>&1')
        shell('echo -e "INFO: Finished creating BWA index for {input.ref}\n"')

def before(value, a):
    # Find first part and return slice before it.
    pos_a = value.find(a)
    if pos_a == -1: return ""
    return value[0:pos_a]

def get_samples_for_bwa(wildcards):
    folder = "datasets/"+wildcards.df+"/reads/"+wildcards.preproc+'/'+wildcards.sample+'/'
    
    # sample name
    sample_name = wildcards.sample
    files = []
    result = ''
    res = []
    
    for fname in os.listdir(folder):
        fl = before(fname, ".fastq")
        last_3 = (fl[len(fl)-3:len(fl)])
        if last_3 == '_R1' or last_3 == '_R2':
            fl = (fl[0:len(fl)-3])
        if sample_name == fl:
            files.append(fname)
            # print(fname)
    i = 0
    if len(files) < 1:
        raise Exception("NO FILES!")
    elif len(files) > 1:
        res.append('')
        res.append('')
        for f in files:
            if '_S.' not in f:
                result += f
                if '_R1.' in f:
                    res[0]=(folder+f)
                elif '_R2.' in f:
                    res[1]=(folder+f)
                if i == 0:
                    result += ' '
                    i += 1
    elif len(files) == 1:
        res.append(folder+str(files[0]))
        result += files[0]
    return res

def get_ref(wildcards):
    wc_str = 'data/ref/index/bwa/{type}/{id_seq_set}/index.sa'
    if wildcards.type == 'assembRaw':
        dfs = ''
        preprocs = ''
        samples = ''
        wc_str = 'processing/assemb/{tool_param}/{dfs}/{preprocs}/{samples}/index/bwa/index.sa'
        df = wildcards.df
        preproc = wildcards.preproc
        
        tool_param, assemb_inp = wildcards.id_seq_set.split('::')
        print(tool_param)
        full_assemb = df+'='+preproc+assemb_inp
        print(full_assemb.split('+'))
        
        spl = full_assemb.split('=')
        dfs = spl[0]
        pr_w_samps = spl[1:] # preprocs with samples 'imp-p136C:p136N'
        for prws in pr_w_samps:
            prws_spl = prws.split('-')
            print(prws_spl)
            preprocs += prws_spl[0]+'='
            samps = prws_spl[1]
            samples += samps+'='
        preprocs=preprocs[0:-1]
        samples=samples [0:-1]
        
        return wc_str.format(tool_param=tool_param, dfs=dfs, preprocs=preprocs, samples=samples)
    else:
        wc_str = nucl_dir+'index/bwa/{type}/{category}/{seq_object}/{seq_set_id}/index.sa'
        return wc_str.format(type=wildcards.type,
                             category=wildcards.category,
                             seq_object=wildcards.seq_object,
                             seq_set_id=wildcards.seq_set_id)
        
        
rule map_on_ref_bwa:
    input:
        r1 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = 'datasets/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        ref_fasta = get_ref
        #'data/ref/index/bwa/{type}/{id_seq_set}/index.sa'
    output:
        sam = 'datasets/{df}/mapped/{preproc}/bwa__{params}/{type}__{category}__{seq_object}/{seq_set_id}/mapped/{sample}.sam'
    log:      'datasets/{df}/mapped/{preproc}/bwa__{params}/{type}__{category}__{seq_object}/{seq_set_id}/mapped/{sample}.log'
    #benchmark: 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/mapped/{sample}.time'
    threads: 6
    run:
        ind_prefix = input.ref_fasta[0:-3]
        shell('({config[bwa.bin]} mem -M -t {threads} {ind_prefix} {input.r1} {input.r2} > {output.sam}) >{log} 2>&1')
        # shell('({config[bwa.bin]} mem -M -t {threads} {ind_prefix} {input.reads} | /srv/common/bin/samtools view -SF 4 -h > {output.sam}) >{log} 2>&1')

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
rule map_on_contig_bwa:
    input:
        reads1 = 'data/reads/{method}/{sample}_R1.trim.nohum.fastq.gz',
        reads2 = 'data/reads/{method}/{sample}_R2.trim.nohum.fastq.gz',
        ref_fasta = 'data/ref/index/bwa/assemb/{contig}-{method}-{tool}/index.sa'
    output:
        sam = 'data/mapped/bwa/assemb/{sample}-{method}_vs_{contig}-{method}-{tool}.sam'
    params: 
        ref_prefix = 'data/ref/index/bwa/assemb/{contig}-{method}-{tool}/index',
        tmp_sam = 'data/mapped/bwa/assemb/{sample}-{method}_vs_{contig}-{method}-{tool}.tmp.sam'
    log: 'log/bwa/map/{sample}-{method}_vs_{contig}-{method}-{tool}.log'
    benchmark: 'time/bwa/map/{sample}-{method}_vs_{contig}-{method}-{tool}.time.txt'
    threads: 20
    run:
        shell('({config[bwa.bin]} mem -M -t {threads} {params.ref_prefix} {input.reads1} {input.reads2} > {params.tmp_sam}) >{log} 2>&1')
        shell('''echo -e "INFO: Adding read groups to SAM file: " \n
               {config[java.bin]} \
                   -jar {config[picard.jar]} AddOrReplaceReadGroups \
                      I={params.tmp_sam} \
                      O={output.sam} \
                      RGID=R_{wildcards.sample} \
                      RGLB=lib1 \
                      RGPL=illumina \
                      RGPU=unit1 \
                      RGSM={wildcards.sample} \n
                      echo -e "\nINFO: Read groups added."''')
        shell('rm {params.tmp_sam}')
        
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