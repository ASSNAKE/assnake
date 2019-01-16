rule init_bam_mapped:
    input:
        sam = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}__{category}__{seq_object}/{id_seq_set}/{what}/{id_sample}.sam'
    output:
        bam = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}__{category}__{seq_object}/{id_seq_set}/{what}/{id_sample}.bam'
    params:
        tmp = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}__{category}__{seq_object}/{id_seq_set}/{what}/{id_sample}.tmp.bam'
    run: 
        shell('echo "{input.sam}"')   
        shell('{config[samtools.bin]} view -bS {input.sam} -o {params.tmp}')
        shell('{config[samtools.bin]} sort {params.tmp} -o {output.bam}')
        shell('{config[samtools.bin]} index {output.bam}')
        shell('rm {params.tmp}')

nucl_dir= config['na_db_dir']
rule coverage:
    input:
        bam = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}/{id_seq_set}/{what}/{id_sample}.bam',
        ref = nucl_dir+"{type}/{fasta}.fa"
    output:
        cov = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}/{id_seq_set}/{what}/{id_sample}.cov'
    run:
        shell('/srv/common/bin/genomeCoverageBed -bga -ibam {input.bam} > {output.cov}')
        
rule coverage_stats:
    input:
        bam = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}__{category}__{seq_object}/{id_seq_set}/{what}/{id_sample}.bam',
        ref = nucl_dir + "{type}/{category}/{seq_object}/{id_seq_set}.fa"
    output:
        stats = 'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}__{category}__{seq_object}/{id_seq_set}/{what}/{id_sample}.bb_stats'
    log:'datasets/{df}/mapped/{preproc}/{mapping_tool}/{type}__{category}__{seq_object}/{id_seq_set}/{what}/{id_sample}_cov_summary.txt'
    run:
        shell('''{config[java.bin]} -ea -Xmx35000m -cp /data6/bio/TFM/soft/bbmap/current/ jgi.CoveragePileup ref={input.ref} in={input.bam} out={output.stats} >{log} 2>&1 \n
        cat {log}''')
             
rule init_bam:
    input:
        sam = 'data/mapped/bwa/assemb/{sample}-{method}_vs_{contig}-{method}-{tool}.sam'
    output:
        bam = 'data/mapped/bwa/assemb/{sample}-{method}_vs_{contig}-{method}-{tool}.bam'
    params:
        tmp = 'data/mapped/bwa/assemb/{sample}-{method}_vs_{contig}-{method}-{tool}.tmp.bam'
    run:    
        shell('{config[samtools.bin]} view -bS {input.sam} -o {params.tmp}')
        shell('{config[samtools.bin]} sort {params.tmp} -o {output.bam}')
        shell('{config[samtools.bin]} index {output.bam}')
