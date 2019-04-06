import os
import pandas as pd

fna_db_dir = '/data5/bio/databases/fna'
contigs_fasta_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__{min_len}.fa')
clean_contigs_fasta_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__{min_len}__no_hum_centr.fa')
contigs_centr_kark_style_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__{min_len}__centr__{params}_krak.tsv')
contigs_centr_class_wc = os.path.join(fna_db_dir, 'assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__{min_len}__centr__{params}_classification.tsv')


contigs_centr_kark_style = contigs_centr_kark_style_wc.format(
    dfs = 'FHM',
    params = 'def',
    preprocs = 'imp__tmtic_def1',
    samples = 'DFM_003_F1_S10:DFM_3F2_S63:D3T3_L1S1_B6_S106',
    min_len = '1000'
)
contigs_centr_class = contigs_centr_class_wc.format(
    dfs = 'FHM',
    params = 'def',
    preprocs = 'imp__tmtic_def1',
    samples = 'DFM_003_F1_S10:DFM_3F2_S63:D3T3_L1S1_B6_S106',
    min_len = '1000'
)
contigs_fasta_loc = contigs_fasta_wc.format(
    dfs = 'FHM',
    params = 'def',
    preprocs = 'imp__tmtic_def1',
    samples = 'DFM_003_F1_S10:DFM_3F2_S63:D3T3_L1S1_B6_S106',
    min_len = '1000'
)
clean_contigs_fasta_loc = clean_contigs_fasta_wc.format(
    dfs = 'FHM',
    params = 'def',
    preprocs = 'imp__tmtic_def1',
    samples = 'DFM_003_F1_S10:DFM_3F2_S63:D3T3_L1S1_B6_S106',
    min_len = '1000'
)

ccc = pd.read_csv(contigs_centr_class, sep = '\t')

human_contigs = list(ccc.loc[ccc['taxID'] == 9606]['readID'])

non_human_contigs = ''

with open(contigs_fasta_loc, 'r') as contigs:
    for i, line in enumerate(contigs):
        if line[1:-1] in human_contigs:
            line = contigs.readline()
        else:
            non_human_contigs += line
            
with open(clean_contigs_fasta_loc, 'w') as clean_contigs:
    clean_contigs.write(non_human_contigs)