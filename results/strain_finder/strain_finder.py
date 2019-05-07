rule sf_preprocess_bam:
    input:
        sam = '{prefix}/{df}/mapped/bwa__{params}/StrainFinder/sf_amphora/{sample}/{preproc}/mapped.sam'
    output:
        sorted_bam = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.bam'
    params:
        filt = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.filt.sam',
        bam = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.presort.bam',
        filt_kumar = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.filt_kumar.sam'
    run:
        shell('echo -e "INFO: Filtering {input.sam}\n"')
        shell('python bin/scripts/StrainFinder/1.filter_sam.py --input {input.sam} > {params.filt}')
        
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
        fa = os.path.join(fna_db_dir, 'StrainFinder/sf_amphora.fa')
    output:
        gf = os.path.join(fna_db_dir, 'StrainFinder/sf_amphora.gf')
    run:
        shell('python bin/scripts/StrainFinder/2.make_gene_file.py --fst {input.fa} --out {output.gf}')
        
rule sf_kpileup:
    input:
        gf = os.path.join(fna_db_dir, 'StrainFinder/sf_amphora.gf'),
        s_bam = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.bam'
    output:
        kpileup = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.kp'
    conda: 'env.yaml'    
    shell: ('perl bin/scripts/StrainFinder/3.kpileup.pl {wildcards.sample} {input.s_bam} {input.gf} 20 0 5 > {output.kpileup}')
        
def get_kp_files(wildcards):
    kp_wc = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample}/mapped.kp'
    sample_list_wc = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample_list}/sample_list.txt'

    sample_list = sample_list_wc.format(df=wildcards.df, 
                           preproc=wildcards.preproc, 
                           prefix=wildcards.prefix, 
                           params=wildcards.params, 
                           sample_list = wildcards.sample_list)
    samples = []
    with open(sample_list, 'r') as sl:
        samples = sl.readlines()

    samples_kp = []
    for s in samples:
        samples_kp.append(kp_wc.format(df=wildcards.df, 
                           preproc=wildcards.preproc, 
                           prefix=wildcards.prefix, 
                           params=wildcards.params,
                           sample=s.strip()))

    return samples_kp    


rule sf_to_numpy:
    input:
        gf = os.path.join(fna_db_dir, 'StrainFinder/sf_amphora.gf'),
        sampe_list =  '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample_list}/sample_list.txt',
        samples = get_kp_files
    output:
        aligments = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample_list}/all_alignments.cPickle'
    conda: '../../envs/python3_std.yaml'    
    shell: ('''which python;\n
        python {config[assnake_install_dir]}/results/strain_finder/scripts/kp2np.py --samples {input.samples} --gene_file {input.gf} --out {output.aligments}''')
        
rule sf_concact:
    input:
        sampe_list =  '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample_list}/sample_list.txt',
        samples = get_kp_files,
        aln = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample_list}/all_alignments.cPickle',
        mapp = os.path.join(fna_db_dir, 'StrainFinder/sf_amphora.mapp'),
    output:
        log = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{sample_list}/all_alignments.log'
    run:
        shell('python {config[assnake_install_dir]}/results/strain_finder/scripts/filter_np.py --aln {input.aln} --map {input.mapp} --samples {input.samples} > {output.log}')

        
rule sf:
    input:
        aln = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{samples}/{org}.np.cPickle'
    output:
        em_out = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{samples}/res/{org}.em',
        otu_out = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{samples}/res/{org}.otu',
        log = '{prefix}/{df}/StrainFinder/bwa__{params}/{preproc}/{samples}/res/{org}.log'
    conda: 'python_env.yaml'
    shell: ('''conda activate /data6/bio/TFM/pipeline/.snakemake/conda/caf26b02;\n
        which python2.7; \n
        /data6/bio/TFM/pipeline/.snakemake/conda/3dfdbd67/bin/python {config[assnake_install_dir]}/results/strain_finder/scripts/StrainFinder.py --aln {input.aln} -N 5 --max_reps 10 --dtol 1 --ntol 2 --max_time 3600 --converge  --em_out {output.em_out} --otu_out {output.otu_out} --log {output.log} --n_keep 3 --force_update --merge_out --msg''')