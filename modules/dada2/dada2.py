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

rule dada2_make_seqtab:
    input: ''
    output: ''
    conda: 'dada2.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/dada2/make_seqtab_wrapper.py') 