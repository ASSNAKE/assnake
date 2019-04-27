CHECKM = config['CHECKM']

rule check_m:
    input:
        # TODO replace with files
        bin_folder = '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/final_contigs__1000__no_hum_centr.fa.metabat-bins'
    output:
        done = '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/checkm.done'
    params:
        wd    = '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/checkm'
    log:           '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/checkm-log.txt'
    benchmark:      '{prefix}/{df}/metabat2/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/checkm-benchmark.txt'
    threads: 20
    conda: 'checkm_env.yaml'
    shell: ('''echo {CHECKM} | checkm data setRoot {CHECKM}; \n
        (checkm lineage_wf -t {threads} -x fa {input.bin_folder} {params.wd}) >{log} 2>&1; \n
        touch {output.done}''')