import yaml

rule dada2_filater_and_trim:
    input: 
        r1 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_R2.fastq.gz',
        params = os.path.join(config['assnake_db'], 'params/dada2/filter_and_trim/{params}.yaml')
    output:
        r1 = '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}_R1.fastq.gz',
        r2 = '{prefix}/{df}/reads/{preproc}__dada2fat_{params}/{sample}/{sample}_R2.fastq.gz'
    run:
        params_str = 'truncLen=c({truncLen_f},{truncLen_r}), maxEE=c({maxEE_f}, {maxEE_r}), truncQ={truncQ}, maxN={maxN}, rm.phix={rm_phix}, compress=TRUE, verbose=TRUE'
        with open(input.params, 'r') as stream:
            try:
                par = yaml.safe_load(stream)
                params_str = params_str.format(**par)
                print(params_str)
            except yaml.YAMLError as exc:
                print(exc)
