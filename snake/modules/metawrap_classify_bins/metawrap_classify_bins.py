install_dir = config['assnake_install_dir']
prune_script = os.path.join(install_dir, 'modules/metawrap_classify_bins/prune_blast_hits.py')
classify_bins = os.path.join(install_dir, 'modules/metawrap_classify_bins/classify_bins.py')
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

rule prune_contigs:
    input: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.megablast_out.raw.fa')
    output: megablast = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.megablast_out.fa'),
        mapping = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.mapping.tax')
    params:os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.megablast_out.pruned.tab')
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


rule taxator_contigs:
    input: 
        megablast = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.megablast_out.fa'),
        mapping = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.mapping.tax')
    output: 
        predictions =os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.predictions.gff3'),
        binned_predictions = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.binned_predictions.txt'),
        contig_taxonomy =  os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.contig_taxonomy.tab'),
    params: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.megablast_out.pruned.tab'),
    conda: '../blastn/env.yaml'
    threads: 1
    shell: ('''export TAXATORTK_TAXONOMY_NCBI={NCBI_TAX}; \n
    cat {input.megablast} | taxator -a megan-lca -t 0.3 -e 0.01 -g {input.mapping} > {output.predictions}; \n
        sort -k1,1 {output.predictions} | binner -n classification -i genus:0.6 > {output.binned_predictions}; \n
            cat {output.binned_predictions} | taxknife -f 2 --mode annotate -s path | grep -v "Could not" | cut -f1,2 > {output.contig_taxonomy};''')

rule extract_predictions_for_bin:
    input: 
        contig_taxonomy =  os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000__no_hum_centr.contig_taxonomy.tab'),
        names = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/bin_by_bin/{binn}/{binn}-contigs.names')
    output: 
        tax = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/bin_by_bin/{binn}/{binn}-contigs.contig_taxonomy.tab')

    # shell: ('''cat {input.contig_taxonomy} | grep "^$line" {input.contig_taxonomy} > {output.tax}''')
    run: 
        names = []
        with open(input.names, 'r') as n:
            names = [nn.strip() for nn in n.readlines()]
        for name in names:
            command_wc = 'grep "^{name}" {input} >> {output} || true'
            command = command_wc.format(name=name, output=output.tax, input=input.contig_taxonomy)
            shell(command)

rule classify_bin:
    input: 
        tax = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/bin_by_bin/{binn}/{binn}-contigs.contig_taxonomy.tab'),
        bin_contigs = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/bin_by_bin/{binn}/{binn}-contigs.fa')
    output:
        os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/bin_by_bin/{binn}/{binn}-bin_taxonomy.tab')
    conda: 'env.yaml'
    threads: 1
    shell: ('''python2.7 {classify_bins} {input.tax} {input.bin_contigs} > {output}''')