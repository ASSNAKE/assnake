rule metaphlan2:
    input:
        r1 = wc_config['fastq_gz_R1_wc'],
        r2 = wc_config['fastq_gz_R2_wc']
    output: 
        o = wc_config['mp2_out']
    params:
        b =  wc_config['mp2_aligment'],
        MPA_PKL = config['MetaPhlAn2']['mpa_v29'],
        BOWTIE2DB = config['MetaPhlAn2']['mpa_v29_bt'],
        INDEX = 'v29_CHOCOPhlAn_201901',
        task_id = config['task_id'] if 'task_id' in config.keys() else None,

    log:  wc_config['mp2_log']
    benchmark:  wc_config['mp2_time']
    threads: 12
    conda: 'env_2.9.12.yaml'
    wrapper: "file://" + os.path.join(config['assnake_install_dir'], 'modules/metaphlan2/wrapper.py')
        
