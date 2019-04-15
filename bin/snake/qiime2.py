rule dada2_qiime2:
    input: demux = '{prefix}/{df}/qiime2/paired-end-demux.qza'
    output: 
        table = '{prefix}/{df}/qiime2/table.qza',
        re_seqs = '{prefix}/{df}/qiime2/rep-seqs.qza',
        denoising_stats = '{prefix}/{df}/qiime2/denoising-statss.qza',
    log:       '{prefix}/{df}/qiime2/dada2_log.txt'
    threads: 20
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
            qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input.demux} \
            --p-trim-left-f 3 \
            --p-trim-left-r 3 \
            --p-trunc-len-f 247 \
            --p-trunc-len-r 247 \
            --o-table {output.table} \
            --o-representative-sequences {output.re_seqs} \
            --o-denoising-stats {output.denoising_stats} \
            --p-n-threads {threads} \
            --verbose >{log} 2>&1''')