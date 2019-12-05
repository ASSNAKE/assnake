rule resanal:
    input:
        sam_fp = 'datasets/{df}/mapped/{preproc}/bwa__def/db__MEGARes__v1.01/megares/mapped/{sample}.sam',
        annot_fp = '/data5/bio/databases/fna/db/MEGARes/v1.01/annotations.csv',
        ref_fp = '/data5/bio/databases/fna/db/MEGARes/v1.01/megares.fa'
    output:
        gene_fp = 'datasets/{df}/resanal/{sample}/{preproc}/gene_fp.tsv',
        group_fp = 'datasets/{df}/resanal/{sample}/{preproc}/group_fp.tsv',
        class_fp = 'datasets/{df}/resanal/{sample}/{preproc}/class_fp.tsv',
        mech_fp = 'datasets/{df}/resanal/{sample}/{preproc}/mech_fp.tsv',
    conda: "/data6/bio/TFM/pipeline/envs/megares.yml"
    shell: ('''resistome \
                   -ref_fp {input.ref_fp} \
                   -sam_fp {input.sam_fp} \
                   -annot_fp {input.annot_fp} \
                   -gene_fp {output.gene_fp} \
                   -group_fp {output.group_fp} \
                   -class_fp {output.class_fp} \
                   -mech_fp {output.mech_fp} \
                   -t 80''')