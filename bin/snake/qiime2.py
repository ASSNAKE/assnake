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

rule summarize_qiime2:
    input: table = '{prefix}/{df}/qiime2/table.qza',
        meta = '{prefix}/{df}/samples_meta.tsv',
    output: 
        table = '{prefix}/{df}/qiime2/table.qzv',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                qiime feature-table summarize \
            --i-table {input.table} \
            --o-visualization {output.table} \
            --m-sample-metadata-file {input.meta}''')

rule tabulate_qiime2:
    input: rep_seqs = '{prefix}/{df}/qiime2/rep-seqs.qza',
    output: 
        rep_seqs = '{prefix}/{df}/qiime2/rep-seqs.qzv',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                qiime feature-table tabulate-seqs \
            --i-data {input.rep_seqs} \
            --o-visualization {output.rep_seqs}''')

rule align_to_tree_mafft_fasttree:
    input: rep_seqs = '{prefix}/{df}/qiime2/rep-seqs.qza',
    output: 
        al = '{prefix}/{df}/qiime2/aligned-rep-seqs.qza',
        masked_al = '{prefix}/{df}/qiime2/masked-aligned-rep-seqs.qza',
        tree = '{prefix}/{df}/qiime2/unrooted-tree.qza',
        rooted_tree = '{prefix}/{df}/qiime2/rooted-tree.qza'

    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                qiime phylogeny align-to-tree-mafft-fasttree \
                --i-sequences {input.rep_seqs} \
                --o-alignment {output.al} \
                --o-masked-alignment {output.masked_al} \
                --o-tree {output.tree} \
                --o-rooted-tree {output.rooted_tree}''')

rule diversity_core:
    input: 
        rooted_tree = '{prefix}/{df}/qiime2/rooted-tree.qza',
        table = '{prefix}/{df}/qiime2/table.qza',
        meta = '{prefix}/{df}/samples_meta.tsv',
    output: 
        done = '{prefix}/{df}/qiime2/core-metrics-results/done'
    params:
        out_folder = '{prefix}/{df}/qiime2/core-metrics-results'

    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                rm -rf {params.out_folder};
                qiime diversity core-metrics-phylogenetic \
                --i-phylogeny {input.rooted_tree} \
                --i-table {input.table} \
                --p-sampling-depth 1109 \
                --m-metadata-file {input.meta} \
                --output-dir {params.out_folder}; \n
                    touch {output.done}''') 

rule diversity_results:
    input: 
        faith_pd_vector = '{prefix}/{df}/qiime2/core-metrics-results/faith_pd_vector.qza',
        evenness_vector = '{prefix}/{df}/qiime2/core-metrics-results/evenness_vector.qza',
        table = '{prefix}/{df}/qiime2/table.qza',
        meta = '{prefix}/{df}/samples_meta.tsv',
    output: 
        faith_pd_group_significance = '{prefix}/{df}/qiime2/faith-pd-group-significance.qzv',
        evenness_group_significance = '{prefix}/{df}/qiime2/evenness-group-significance.qzv',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                qiime diversity alpha-group-significance \
                --i-alpha-diversity {input.faith_pd_vector} \
                --m-metadata-file {input.meta} \
                --o-visualization {output.faith_pd_group_significance};\n

                qiime diversity alpha-group-significance \
                --i-alpha-diversity {input.evenness_vector} \
                --m-metadata-file {input.meta} \
                --o-visualization {output.evenness_group_significance}''') 

QIIME2_SILVA99 = config['QIIME2_SILVA99']

rule classify_qiime2:
    input: 
        reads          = '{prefix}/{df}/qiime2/rep-seqs.qza',
    output: 
        classification = '{prefix}/{df}/qiime2/taxonomy.qza',
        # class_v        = '{prefix}/{df}/qiime2/taxonomy.qzv',
    threads: 6
    log: '{prefix}/{df}/qiime2/taxonomy.log',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                (qiime feature-classifier classify-sklearn \
                    --i-classifier {QIIME2_SILVA99} \
                    --i-reads {input.reads} \
                    --o-classification {output.classification} \
                    --verbose \
                    --p-n-jobs {threads}) > {log} 2>&1;\n
                ''') 

                # qiime metadata tabulate \
                #     --m-input-file {output.classification} \
                #     --o-visualization {output.class_v}
QIIME2_GG99 = config['QIIME2_GG99']

rule classify_qiime2_gg:
    input: 
        reads          = '{prefix}/{df}/qiime2/rep-seqs.qza',
    output: 
        classification = '{prefix}/{df}/qiime2/gg/taxonomy.qza',
        # class_v        = '{prefix}/{df}/qiime2/taxonomy.qzv',
    threads: 6
    log: '{prefix}/{df}/qiime2/gg/taxonomy.log',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                (qiime feature-classifier classify-sklearn \
                    --i-classifier {QIIME2_GG99} \
                    --i-reads {input.reads} \
                    --reads-per-batch 1000 \
                    --o-classification {output.classification} \
                    --verbose \
                    --p-n-jobs {threads}) > {log} 2>&1;\n
                ''') 

rule classify_qiime2_gg_meta:
    input: 
        classification = '{prefix}/{df}/qiime2/gg/taxonomy.qza',
    output: 
        class_v        = '{prefix}/{df}/qiime2/gg/taxonomy.qzv',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                qiime metadata tabulate \
                     --m-input-file {input.classification} \
                     --o-visualization {output.class_v}
                ''') 

rule taxa_bar:
    input: 
        classification = '{prefix}/{df}/qiime2/gg/taxonomy.qza',
        table = '{prefix}/{df}/qiime2/table.qza',
        meta = '{prefix}/{df}/samples_meta.tsv',
    output: 
        bar        = '{prefix}/{df}/qiime2/gg/taxa-bar-plots.qzv',
    shell: ('''source /data4/bio/fedorov/miniconda3/bin/activate qiime2-2019.1; \n
                    qiime taxa barplot \
                --i-table {input.table} \
                --i-taxonomy {input.classification} \
                --m-metadata-file {input.meta} \
                --o-visualization {output.bar}''')