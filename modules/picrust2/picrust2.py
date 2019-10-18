
fasta_from_seqtab_script = os.path.join(config['assnake_install_dir'], 'modules/picrust2/fasta_from_seqtab.R')
rule fasta_from_seqtab:
    input: '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.rds'
    output: '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim_seqs.fa'
    conda: '../dada2/dada2.yaml'
    shell: "Rscript {fasta_from_seqtab_script} '{input}' '{output}'"

db = config['picrust2_pro']
rule picrust2_tree:
    input: '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim_seqs.fa'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/placed_seqs.tre'
    params: '{prefix}/{df}/dada2/{sample_set}/picrust2/placement_working'
    conda: 'picrust2_env.yaml'
    threads: 20
    shell: 'place_seqs.py -s {input} -o {output} -p {threads} --ref_dir {db} --intermediate {params}'

rule picrust2_16s:
    input: '{prefix}/{df}/dada2/{sample_set}/picrust2/placed_seqs.tre'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/marker_nsti_predicted.tsv.gz'
    conda: 'picrust2_env.yaml'
    threads: 20
    shell: 'hsp.py -i 16S -t {input} -o {output} -p {threads} -n'

rule picrust2_ec:
    input: '{prefix}/{df}/dada2/{sample_set}/picrust2/placed_seqs.tre'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_predicted.tsv.gz'
    conda: 'picrust2_env.yaml'
    threads: 20
    shell: 'hsp.py -i EC -t {input} -o {output} -p {threads}'

rule picrust2_ko:
    input: '{prefix}/{df}/dada2/{sample_set}/picrust2/placed_seqs.tre'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/KO_predicted.tsv.gz'
    conda: 'picrust2_env.yaml'
    threads: 20
    shell: 'hsp.py -i KO -t {input} -o {output} -p {threads}'

rule metagenome_pipeline_ec:
    input: 
        ec = '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_predicted.tsv.gz',
        ma = '{prefix}/{df}/dada2/{sample_set}/picrust2/marker_nsti_predicted.tsv.gz',
        co = '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.tsv'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/mg_ec.done'
    params: '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_metagenome_out'
    conda: 'picrust2_env.yaml'
    threads: 20
    shell: 'metagenome_pipeline.py -i {input.co} \
                       -m {input.ma} \
                       -f {input.ec} \
                       -o {params}; touch {output}'
    
rule metagenome_pipeline_ko:
    input: 
        ec = '{prefix}/{df}/dada2/{sample_set}/picrust2/KO_predicted.tsv.gz',
        ma = '{prefix}/{df}/dada2/{sample_set}/picrust2/marker_nsti_predicted.tsv.gz',
        co = '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.tsv'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/mg_ko.done'
    params: '{prefix}/{df}/dada2/{sample_set}/picrust2/KO_metagenome_out'
    conda: 'picrust2_env.yaml'
    shell: 'metagenome_pipeline.py -i {input.co} \
                       -m {input.ma} \
                       -f {input.ec} \
                       -o {params}; touch {output}'

rule metagenome_pipeline_ec_paramed:
    input: 
        ec = '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_predicted.tsv.gz',
        ma = '{prefix}/{df}/dada2/{sample_set}/picrust2/marker_nsti_predicted.tsv.gz',
        co = '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.tsv'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/mg_ec__david1.done'
    params: '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_metagenome_out__david1'
    conda: 'picrust2_env.yaml'
    threads: 20
    shell: 'metagenome_pipeline.py -i {input.co} \
                        -m {input.ma} \
                        -f {input.ec} \
                        --min_reads 10 \
                        --min_samples 18 \
                        -o {params}; touch {output}'
rule metagenome_pipeline_ko_paramed:
    input: 
        ec = '{prefix}/{df}/dada2/{sample_set}/picrust2/KO_predicted.tsv.gz',
        ma = '{prefix}/{df}/dada2/{sample_set}/picrust2/marker_nsti_predicted.tsv.gz',
        co = '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.tsv'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/mg_ko__david1.done'
    params: '{prefix}/{df}/dada2/{sample_set}/picrust2/KO_metagenome_out__david1'
    conda: 'picrust2_env.yaml'
    shell: 'metagenome_pipeline.py -i {input.co} \
                        -m {input.ma} \
                        -f {input.ec} \
                        --min_reads 10 \
                        --min_samples 18 \
                        -o {params}; touch {output}'

rule pathway_pipeline_paramed:
    input: ec = '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_metagenome_out__david1/pred_metagenome_unstrat.tsv.gz'
    output: done = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways__david1.done',
        a = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out__david1/path_abun_unstrat.tsv.gz',
    params: pw_out = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out__david1',
        pw_work = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_working__david1'
    conda:'picrust2_env.yaml'
    threads: 20
    shell: 'pathway_pipeline.py -i {input.ec} \
                    -o {params.pw_out} \
                    --intermediate {params.pw_work} \
                    -p {threads}; touch {output.done}'

rule add_desc_pathway_paramed:
    input:  
        a = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out__david1/path_abun_unstrat.tsv.gz',
        d = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways__david1.done'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out__david1/path_abun_unstrat_descrip.tsv.gz'
    conda:'picrust2_env.yaml'
    shell: 'add_descriptions.py -i {input.a} -m METACYC \
                    -o {output}'


rule pathway_pipeline:
    input: ec = '{prefix}/{df}/dada2/{sample_set}/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz'
    output: done = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways.done'
    params: pw_out = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out',
        pw_work = '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_working'
    conda:'picrust2_env.yaml'
    threads: 20
    shell: 'pathway_pipeline.py -i {input.ec} \
                    -o {params.pw_out} \
                    --intermediate {params.pw_work} \
                    -p {threads}; touch {output.done}'
                    

rule add_desc_pathway:
    input:  '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out/path_abun_unstrat.tsv.gz'
    output: '{prefix}/{df}/dada2/{sample_set}/picrust2/pathways_out/path_abun_unstrat_descrip.tsv.gz'
    conda:'picrust2_env.yaml'
    shell: 'add_descriptions.py -i {input} -m METACYC \
                    -o {output}'