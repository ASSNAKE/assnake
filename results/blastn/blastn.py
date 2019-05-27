BLASTDB = config['BLAST_NCBI_NT']

rule classify_bins:
    input: '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.fa'
    output:'{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.megablast_out.raw.tab'
    conda: 'env.yaml'
    threads: 10
    shell: ('''blastn -task megablast -num_threads {threads} \
            -db {BLASTDB}/nt \
            -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length' \
            -query {input} > {output}''')

rule classify_contigs:
    input: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.fa')
    output: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.megablast_out.raw.fa')
    conda: 'env.yaml'
    threads: 12
    shell: ('''blastn -task megablast -num_threads {threads} \
            -db {BLASTDB}/nt \
            -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length' \
            -query {input} > {output}''')

rule classify_contigs_2:
    input: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000.fa')
    output: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000.megablast_out.raw.fa')
    conda: 'env.yaml'
    threads: 12
    shell: ('''blastn -task megablast -num_threads {threads} \
            -db {BLASTDB}/nt \
            -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length' \
            -query {input} > {output}''')