import yaml

rule dada2_filter_and_trim:
    input: 
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        params = os.path.join(config['assnake_db'], 'params/dada2/filter_and_trim/{params}.yaml')
    output:
        r1 = '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}_R2.fastq.gz'
    log: '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}.log'
    conda: 'dada2.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/dada2/filter_trim_wrapper.py')

rule dada2_learn_errors:
    input: 
        samples_list = os.path.join(config['dada2_dir'], '{sample_set}', 'samples.tsv'),
        params = os.path.join(config['assnake_db'], 'params/dada2/learn_errors/{params}.yaml')
    output:
        err          = os.path.join(config['dada2_dir'], '{sample_set}/{params}/err{strand}.rds')
    log:               os.path.join(config['dada2_dir'], '{sample_set}/{params}/err{strand}.log')
    conda: 'dada2.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/dada2/learn_errors_wrapper.py')


derep_dada_merge_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/derep_dada_merge.R')
rule dada2_derep_dada_merge:
    input: 
        r1     = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2     = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        errF   = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/errR1.rds'),
        errR   = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/errR2.rds'),
        params = os.path.join(config['assnake_db'], 'params/dada2/merge/{params}.yaml')
    output:
        merged = '{prefix}/{df}/reads/{preproc}/{sample}/dada2/{sample_set}/{err_params}/merged_{params}.rds',
        stats  = '{prefix}/{df}/reads/{preproc}/{sample}/dada2/{sample_set}/{err_params}/merged_{params}.stats',
    log: '{prefix}/{df}/reads/{preproc}/{sample}/dada2/{sample_set}/{err_params}/log_{params}.txt'
    conda: 'dada2.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/dada2/derep_merge_wrapper.py') 

rule dada2_derep_infer_pooled:
    input: 
        samples_list = os.path.join(config['dada2_dir'], '{sample_set}', 'samples.tsv'),
        err   = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/err{strand}.rds'),
    output:
        infered = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/dada{strand}.rds'),
        derep = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/derep{strand}.rds'),
    conda: 'dada2.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/dada2/infer_pooled_wrapper.py') 

merge_pooled_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/merge_pooled.R')
rule dada2_merge_pooled:
    input: 
        dada_1 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/dadaR1.rds'),
        dada_2 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/dadaR2.rds'),
        derep_1 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/derepR1.rds'),
        derep_2 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/derepR2.rds'),
    output:
        mergers = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/mergers_{len}.rds'),
    wildcard_constraints:    
        sample_set="[\w\d_-]+",
        err_params="[\w\d_-]+"
    conda: 'dada2.yaml'
    shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript {merge_pooled_script} '{input.derep_1}' '{input.derep_2}'  '{input.dada_1}' '{input.dada_2}' '{output.mergers}';''') 


rule dada2_derep_infer_pooled_sub:
    input: 
        samples_list = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}', 'samples.tsv'),
        err   = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/err{strand}.rds'),
    output:
        infered = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/dada{strand}.rds'),
        derep = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/derep{strand}.rds'),
    wildcard_constraints:    
        err_params="[\w\d_-]+"
    conda: 'dada2.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/dada2/infer_pooled_wrapper.py') 

rule dada2_merge_pooled_sub:
    input: 
        dada_1 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/dadaR1.rds'),
        dada_2 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/dadaR2.rds'),
        derep_1 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/derepR1.rds'),
        derep_2 = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/derepR2.rds'),
    output:
        mergers = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/mergers_{len}.rds'),
    conda: 'dada2.yaml'
    shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript {merge_pooled_script} '{input.derep_1}' '{input.derep_2}'  '{input.dada_1}' '{input.dada_2}' '{output.mergers}';''') 

derep_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/derep.R')
# rule dada2_derep:
#     input:
#         samples_list = os.path.join(config['dada2_dir'], '{sample_set}/', 'samples.tsv'),
#         err   = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/err{strand}.rds'),
#     output:
#         derep = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/derep{strand}.rds'),
#     conda: 'dada2.yaml'
#     shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
#         Rscript {derep_script} '{input.samples_list}' '{wildcards.strand}' '{output.derep}';''')

make_seqtab_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/make_seqtab.R')
rule dada2_make_seqtab:
    input: mergers = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/mergers_{len}.rds'),
    output: os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/seqtab_{len}.rds')
    conda: 'dada2.yaml'
    wildcard_constraints:    
        sample_set="[\w\d_-]+",
        err_params="[\w\d_-]+"
    shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript {make_seqtab_script} '{input.mergers}' '{output}';''') 

seqtab_nochim_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/seqtab_nochim.R')
rule dada2_nochim:
    input: '{prefix}/{df}/dada2/{sample_set}/seqtab.rds'
    output: '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.rds'
    conda: 'dada2.yaml'
    wildcard_constraints:    
        sample_set="[\w\d_-]+",
        err_params="[\w\d_-]+"
    shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript {seqtab_nochim_script} '{input}' '{output}';''') 

rule dada2_make_seqtab_sub:
    input: mergers = os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/mergers_{len}.rds'),
    output: os.path.join(config['dada2_dir'], '{sample_set}/{err_params}/{sub}/seqtab_{len}.rds')
    conda: 'dada2.yaml'
    shell: ('''export LANG=en_US.UTF-8;\nexport LC_ALL=en_US.UTF-8;\n
        Rscript {make_seqtab_script} '{input.mergers}' '{output}';''') 

assign_taxa_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/assign_taxa.R')
rule dada2_assign_taxa:
    input: '{prefix}/{df}/dada2/{sample_set}/seqtab_nochim.rds'
    output: '{prefix}/{df}/dada2/{sample_set}/taxa.rds'
    conda: 'dada2.yaml'
    shell: ("Rscript {assign_taxa_script} '{input}' '{output}'")

make_tree_script = os.path.join(config['assnake_install_dir'], 'modules/dada2/scripts/make_tree.R')
rule make_tree:
    input: '{prefix}/{df}/dada2/{sample_set}/seqtab_{mod}.rds'
    output: tree = '{prefix}/{df}/dada2/{sample_set}/tree_{mod}.rds',
        al = '{prefix}/{df}/dada2/{sample_set}/aligment_{mod}.rds'
    conda: 'for_tree.yaml'
    shell: ("Rscript {make_tree_script} '{input}' '{output.tree}' '{output.al}'")

