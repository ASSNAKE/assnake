rule humann2:
    input:
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        mp2 = '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/{sample}.mp2'
    output:
        gf = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_genefamilies.tsv',
        pc = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_pathcoverage.tsv',
        pa = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_pathabundance.tsv'
    params:
        wd =     '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/',
        merged = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}.fastq.gz'
    threads: 48
    conda: "../../envs/humann2.yml"
    shell: ('''cat {input.r1} {input.r2} > {params.merged} \n 
               humann2 --taxonomic-profile {input.mp2} \
               --protein-database /data5/bio/databases/humann2/uniref90/uniref \
               --nucleotide-database /data5/bio/databases/humann2/chocophlan/chocophlan \
               --input {params.merged} --output {params.wd} --threads {threads} \n
               rm {params.merged}''') 

mapping = '/data5/bio/runs-jeniaole/tools/humann/data/humann2/kegg/kegg_idmapping.tsv'
pathway = '/data5/bio/runs-jeniaole/tools/humann/data/humann2/kegg/keggc'
custom_db = '/data5/bio/databases/humann2/custom/KEGG_FHM/kegg_fhm_db'
rule humann2_custom_kegg:
    input:
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        mp2 = '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/{sample}.mp2'
    output:
        gf = '{prefix}/{df}/humann2/KEGG/{sample}/{preproc}/{sample}_genefamilies.tsv',
        pc = '{prefix}/{df}/humann2/KEGG/{sample}/{preproc}/{sample}_pathcoverage.tsv',
        pa = '{prefix}/{df}/humann2/KEGG/{sample}/{preproc}/{sample}_pathabundance.tsv'
    params:
        wd =     '{prefix}/{df}/humann2/KEGG/{sample}/{preproc}/',
        merged = '{prefix}/{df}/humann2/KEGG/{sample}/{preproc}/{sample}.fastq.gz'
    threads: 20
    conda: "../../envs/humann2.yml"
    shell: ('''cat {input.r1} {input.r2} > {params.merged} \n 
               humann2 --taxonomic-profile {input.mp2} \
               --id-mapping {mapping} --pathways-database {pathway} \
               --protein-database {custom_db} --bypass-nucleotide-search \
               --input {params.merged} --output {params.wd} --threads {threads} \n
               rm {params.merged}''') 
        
rule humann2_split_stratified_table:
    input:
        pa = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_pathabundance_norm.tsv'
    output:
        strat = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_pathabundance_norm_stratified.tsv',
        unstrat = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_pathabundance_norm_unstratified.tsv',
    params:
        wd =     '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/',
    conda: "../../envs/humann2.yml"
    shell: ('''humann2_split_stratified_table -i {input.pa} -o {params.wd}''')
        
        
rule humann2_regroup:
    input:
        gf = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_genefamilies.tsv',
    output:
        gf = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}__{groups}.tsv',
    conda: "../../envs/humann2.yml"
    shell: ('''humann2_regroup_table -i {input.gf} -o {output} --custom /data5/bio/databases/humann2/ut_mapping/utility_mapping/{wildcards.groups}.txt.gz''') 
        
rule humann2_norm:
    input:
        gf = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_{groups}.tsv',
    output:
        norm = '{prefix}/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}_{groups}_norm.tsv'
    conda: "../../envs/humann2.yml"
    shell: ('''humann2_renorm_table --input {input.gf} --units relab --output {output.norm}''') 

        
# /data5/bio/runs-jeniaole/tools/humann/data/humann2/kegg
# /data5/bio/runs-jeniaole/tools/kegg/ftp.bioinformatics.jp/kegg
# 1 - Join all of the taxonomic profiles
# humann2_join_tables -i ./ --file_name .mp2 -s -o joined
# 2 - Reduce this file into a taxonomic profile that represents the maximum abundances from all of the samples in your set
# humann2_reduce_table --input joined_taxonomic_profile.tsv --output max_taxonomic_profile.tsv --function max --sort-by level
# 3 - Create a custom Kegg database for your data set, with genus-specific taxonomic limitation, using your joint taxonomic profile
# humann2_build_custom_database --input /data5/bio/runs-jeniaole/tools/humann/data/humann2/kegg/genes.pep --output kegg_fhm_db --id-mapping /data5/bio/runs-jeniaole/tools/humann/data/humann2/kegg/kegg_idmapping.tsv --format diamond --taxonomic-profile max_taxonomic_profile.tsv
# humann2 --input $SAMPLE.fastq --output $OUTPUT_DIR --id-mapping legacy_kegg_idmapping.tsv --pathways-database humann1/data/keggc --protein-database custom_database --bypass-nucleotide-search