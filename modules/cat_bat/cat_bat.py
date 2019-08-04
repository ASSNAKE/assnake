cat_bat_db = config['CAT_BAT_DB']
cat_bat_taxa = config['CAT_BAT_TAXA']

rule cat_bat_CAT:
    input:
        exported = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/all_bins.done')
    output:
        done = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/cat_bins.done')
    params:
        bin_folder = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/all_bins/'),
        wd = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/all_bins_cat/'),
        prefix = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/all_bins_cat/out.BAT'),
    threads: 20
    conda: 'cat_bat_env.yaml'
    shell: ('''mkdir -p {params.wd}; \n
        CAT bins -b {params.bin_folder} -d {cat_bat_db} -t {cat_bat_taxa} --bin_suffix .fa --out_prefix {params.prefix} -n {threads}; \n
            touch {output.done}''')

rule cat_bat_contigs:
    input:
        fa     = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000.fa')
    output:
        done   = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000_CAT/cat.done'),
        classification = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000_CAT/out.CAT.contig2classification.txt')
    params:
        wd     = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000_CAT/'),
        prefix = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000_CAT/out.CAT'),
    threads: 24
    conda: 'cat_bat_env.yaml'
    shell: ('''mkdir -p {params.wd}; \n
            CAT contigs -c {input.fa} -d {cat_bat_db} -t {cat_bat_taxa} --out_prefix {params.prefix} -n {threads}; \n
            touch {output.done}''')

rule cat_bat_contigs_new:
    input:
        fa     = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{sample_set}/final_contigs__1000.fa')
    output:
        done   = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{sample_set}/final_contigs__1000_CAT/cat.done'),
        classification = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{sample_set}/final_contigs__1000_CAT/out.CAT.contig2classification.txt')
    params:
        wd     = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{sample_set}/final_contigs__1000_CAT/'),
        prefix = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{sample_set}/final_contigs__1000_CAT/out.CAT'),
    threads: 24
    conda: 'cat_bat_env.yaml'
    shell: ('''mkdir -p {params.wd}; \n
            CAT contigs -c {input.fa} -d {cat_bat_db} -t {cat_bat_taxa} --out_prefix {params.prefix} -n {threads}; \n
            touch {output.done}''')

rule add_names2:
    input:
        classification = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000_CAT/out.CAT.contig2classification.txt')
    output:
        names = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__1000_CAT/contig2classification_names_official.txt')
    conda: 'cat_bat_env.yaml'
    shell: ('''CAT add_names -i {input.classification} -o {output.names} -t {cat_bat_taxa} --only_official''')

rule add_names:
    input:
        done           = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/cat_bins.done'),
        classification = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/all_bins_cat/out.BAT.bin2classification.txt')
    output:
        names = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/{collection}/all_bins_cat/bin2classification_names_official.txt')
    threads: 20
    conda: 'cat_bat_env.yaml'
    shell: ('''CAT add_names -i {input.classification} -o {output.names} -t {cat_bat_taxa} --only_official''')


    