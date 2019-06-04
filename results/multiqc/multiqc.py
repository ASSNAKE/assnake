rule multiqc_fastqc:
    input:
        sample_list = '{prefix}/{df}/multiqc/{sample_set}/{strand}/samples.list'
    output:
        multiqc_report = '{prefix}/{df}/multiqc/{sample_set}/{strand}/multiqc_report.html'
    params:
        wd = '{prefix}/{df}/multiqc/{sample_set}/{strand}/'
    conda: "multiqc.yaml"
    shell: ("export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8; multiqc --file-list {input.sample_list} -o {params.wd}")