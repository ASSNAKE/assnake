def get_fastqc_files_for_multiqc(wildcards):
    fastqc_files = []
    with open('{prefix}/{df}/multiqc/{sample_set}/{strand}/samples.list'.format(**wildcards), 'r') as sample_list:
        fastqc_files = [s.strip() for s in sample_list.readlines()]
    return fastqc_files

rule multiqc_fastqc:
    input:
        sample_list = '{prefix}/{df}/multiqc/{sample_set}/{strand}/samples.list',
        # samples = get_fastqc_files_for_multiqc
    output:
        multiqc_report = '{prefix}/{df}/multiqc/{sample_set}/{strand}/multiqc_report.html'
    params:
        wd = '{prefix}/{df}/multiqc/{sample_set}/{strand}/'
    conda: "multiqc.yaml"
    shell: ("export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8; multiqc --file-list {input.sample_list} -o {params.wd}")