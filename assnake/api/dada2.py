import pandas as pd
import os
from assnake.api.dataset import Dataset

def get_full_track(fs_prefix, df, preproc):
    '''
    Loads infromation about reads number on all steps of dada2 analysis.
    Tracking script from R should be runned first.
    '''
    tr = Dataset(df)
    
#     preproc = 'raw__dada2fat_jgun_b2'
    raw = pd.DataFrame(tr.sample_sets['raw'])
    raw.index = raw['fs_name']

    track_f = pd.DataFrame(tr.sample_sets[preproc])
    track_f = track_f.rename({'reads': "trimmed"}, axis = 1)
    track_f.index = track_f['fs_name']
    track_f['raw'] = raw['reads']
    track_f = track_f.drop(['df', 'prefix', 'preproc', 'fs_name'], axis=1)
    track_f = track_f//2
    
    track_loc = os.path.join(fs_prefix,df,'dada2/pooled/track.tsv')
    track = pd.read_csv(track_loc, sep= '\t')

    track.index = track['Unnamed: 0'].str.replace('_R1.fastq.gz', '')
    track.index.name = 'fs_name'
    track = track.drop(['Unnamed: 0'], axis=1)
    
    mm = track_f.merge(track, left_index = True, right_index = True)
    return mm

def rename_seqs_to_asvs(otu_table, otu_table_renamed_loc, asvs_fasta_loc):
    ''' 
    Because by default from dada2 we have indexes as full sequences we need to give them new human-readable names,
    and save this mapping to fasta file.
    '''
    asv_sequences = list(otu_table.index)
    seq_to_asv_name = {}
    for i, asv_seq in enumerate(asv_sequences):
        seq_to_asv_name.update({asv_seq: 'ASV'+str(i)})
    otu_table_reind = otu_table.rename(index=seq_to_asv_name)
        
    # otu_table_reind_loc = os.path.join(fs_prefix, df, 'dada2/pooled/seqtab_nochim_renamed.tsv')
    otu_table_reind.to_csv(otu_table_renamed_loc, sep='\t')

    asvs_str = ''
    for k, v in seq_to_asv_name.items():
        asvs_str = asvs_str + '>' + v + "\n" + k + '\n'
    asvs_str = asvs_str[0:-1] #remove last \n

    # asvs_fa_loc = os.path.join(fs_prefix, df, 'dada2/pooled/asvs.fa')
    with open(asvs_fasta_loc, 'x') as asvs_fa_file:
        asvs_fa_file.write(asvs_str)
