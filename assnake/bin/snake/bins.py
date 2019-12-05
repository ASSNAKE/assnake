def get_bins(wildcards):
    bin_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/bin_by_bin')

rule exported_bins_to_folder:
    input:  os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/bins_summary.txt'),
    output: os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/all_bins.done')
    params: 
        all_bins_f = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/all_bins/'), 
        bins_f     = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{collection}/bin_by_bin/*/*-contigs.fa')
    run:
        bins = glob(params.bins_f)
        os.mkdir(params.all_bins_f)
        for b in bins:
            b_name = b.split('/')[-2]
            new_loc = os.path.join(params.all_bins_f, b_name + '.fa')
            shutil.copyfile(b, new_loc)
        shell('''touch {output}''')