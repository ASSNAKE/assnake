rule sf_preprocess_bam:
    input:
        sam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.sam'
    output:
        sorted_bam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.sam'
    params:
        filt = 'data/snv/StrainFinder/{type}/{contig}/{sample}-{method}.filt.sam',
        bam = 'data/snv/StrainFinder/{type}/{contig}/{sample}-{method}.presort.bam',
        filt_kumar = 'data/snv/StrainFinder/{type}/{contig}/{sample}-{method}.filt_kumar.sam'
    run:
        shell('echo -e "INFO: Filtering {input.sam}\n"')
        shell('{config[python.bin]} bin/scripts/StrainFinder/1.filter_sam.py --input {input.sam} > {params.filt}')
        
        shell('echo -e "INFO: Converting {params.filt} to BAM\n"')
        shell('{config[samtools.bin]} view -bS -F 4 -o {params.bam} {params.filt}')
        
        shell('echo -e "INFO: Sorting {params.bam}\n"')
        shell('{config[samtools.bin]} sort {params.bam} -o {output.sorted_bam}')
        
        shell('echo -e "INFO: Indexing {output.sorted_bam}\n"')
        shell('{config[samtools.bin]} index {output.sorted_bam}')
        
        shell('echo -e "INFO: Cleanup\n"')
        shell('rm {params.bam} {params.filt}')

rule sf_make_gf:
    input:
        fa = "data/ref/{type}/{contig}.fa"
    output:
        gf = 'data/ref/{type}/{contig}.sfgf'
    run:
        shell('{config[python.bin]} bin/scripts/StrainFinder/2.make_gene_file.py --fst {input.fa} --out {output.gf}')
        
rule sf_kpileup:
    input:
        gf = 'data/ref/{type}/{id_seq_set}.sfgf',
        s_bam = 'datasets/{df}/mapped/{preproc}/bwa/{type}/{id_seq_set}/filtered/{sample}.bam'
    output:
        kpileup = 'datasets/{df}/snv/{preproc}/StrainFinder/{type}/{id_seq_set}/{sample}.kp'
    run:
        shell('perl bin/scripts/StrainFinder/3.kpileup.pl {wildcards.sample} {input.s_bam} {input.gf} 20 0 5 > {output.kpileup}')
        
def get_kp_files(wildcards):
    prefix = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/'
    prefix = prefix.format(df=wildcards.df, 
                           preproc=wildcards.preproc, 
                           seq_type=wildcards.seq_type, 
                           id_seq_set=wildcards.id_seq_set)
    
    samples = wildcards.samples.split(':')
    samples_kp = []
    for s in samples:
        samples_kp.append(prefix + s + '.kp')
    return samples_kp    


rule sf_to_numpy:
    input:
        gf = 'data/ref/{seq_type}/{id_seq_set}.sfgf',
        samples = get_kp_files
    output:
        aligments = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/all_alignments.cPickle'
    run:
        shell('{config[python.bin]} bin/scripts/StrainFinder/4.kp2np.py --samples {input.samples} --gene_file {input.gf} --out {output.aligments}')
        
rule sf_concact:
    input:
        samples = get_kp_files,
        aln = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/all_alignments.cPickle',
        mapp = 'data/ref/{seq_type}/{id_seq_set}.map'
    output:
        log = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/all_alignments.log'
    run:
        shell('{config[python.bin]} bin/scripts/StrainFinder/5.filter_np.py --aln {input.aln} --map {input.mapp} --samples {input.samples} > {output.log}')

        
rule sf:
    input:
        aln = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/BACT_546.np.cPickle'
    output:
        em_out = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/res/BACT_546.em',
        otu_out = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/res/BACT_546.otu',
        log = 'datasets/{df}/snv/{preproc}/StrainFinder/{seq_type}/{id_seq_set}/{samples}/res/BACT_546.log'
    run:
        shell(' /home/fedorov/miniconda3/envs/sf/bin/python bin/scripts/StrainFinder/StrainFinder.py --aln {input.aln} -N 5 --max_reps 10 --dtol 1 --ntol 2 --max_time 3600 --converge  --em_out {output.em_out} --otu_out {output.otu_out} --log {output.log} --n_keep 3 --force_update --merge_out --msg')