install_dir = config['assnake_install_dir']
prune_script = os.path.join(install_dir, 'results/metawrap_classify_bins/prune_blast_hits.py')
classify_bins = os.path.join(install_dir, 'results/metawrap_classify_bins/classify_bins.py')
NCBI_TAX =config['NCBI_TAX']
NCBI_NODES = os.path.join(config['NCBI_TAX'],'nodes.dmp')

rule prune:
    input: '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.megablast_out.raw.tab'
    output: megablast = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.megablast_out.tab',
        mapping = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.mapping.tax'
    params:'{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.megablast_out.pruned.tab'
    conda: 'env.yaml'
    threads: 1
    shell: ('''python2.7 {prune_script} {NCBI_NODES} {input} > {params}; \n
                cat {params} | cut -f1,2,3,4,5,7,8,9,10,11,12 > {output.megablast}; \n
                cat {params} | cut -f5,6 > {output.mapping}''')

rule taxator:
    input: 
        megablast = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.megablast_out.tab',
        mapping = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.mapping.tax'
    output: predictions = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.predictions.gff3',
        binned_predictions = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.binned_predictions.txt',
        contig_taxonomy =  '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.contig_taxonomy.tab',
    params:'{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.megablast_out.pruned.tab'
    conda: '../blastn/env.yaml'
    threads: 1
    shell: ('''export TAXATORTK_TAXONOMY_NCBI={NCBI_TAX}; \n
    cat {input.megablast} | taxator -a megan-lca -t 0.3 -e 0.01 -g {input.mapping} > {output.predictions}; \n
        sort -k1,1 {output.predictions} | binner -n classification -i genus:0.6 > {output.binned_predictions}; \n
            cat {output.binned_predictions} | taxknife -f 2 --mode annotate -s path | grep -v "Could not" | cut -f1,2 > {output.contig_taxonomy};''')

rule classify_bin:
    input: 
        tax = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.contig_taxonomy.tab',
        bin_contigs = '{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-contigs.fa'
    output:'{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/SUMMARY/bin_by_bin/{bin}/{bin}-bin_taxonomy.tab'
    conda: 'env.yaml'
    threads: 1
    shell: ('''python2.7 {classify_bins} {input.tax} {input.bin_contigs} > {output}''')