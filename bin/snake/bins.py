def get_bins(wildcards):
    bin_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/conocot_anvio5_def/bin_by_bin')

rule exported_bins_to_folder:
    input: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/conocot_anvio5_def/export.done'),
    output: os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/conocot_anvio5_def/all_bins.done')
    params: 
        all_bins_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/conocot_anvio5_def/all_bins/'), 
        bins_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/conocot_anvio5_def/bin_by_bin/*/*-contigs.fa')
    run:
        bins = glob(params.bins_f)
        os.mkdir(params.all_bins_f)
        for b in bins:
            b_name = b.split('/')[-2]
            new_loc = os.path.join(params.all_bins_f, b_name + '.fa')
            shutil.copyfile(b, new_loc)
        shell('''touch {output}''')