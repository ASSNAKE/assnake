import os
import pandas as pd


# contigs_centr_class
# contigs_fasta_loc
# clean_contigs_fasta_loc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--classification', type=str, default='',
                        help='Centrifuge classification')
    parser.add_argument('--contigs', type=str, default='',
                        help='Fasta with contigs')
    parser.add_argument('--clean', type=str, default='',
                        help='File to write clean contigs')
    args = parser.parse_args()

    ccc = pd.read_csv(args.classification, sep = '\t')

    human_contigs = list(ccc.loc[ccc['taxID'] == 9606]['readID'])

    non_human_contigs = ''

    with open(args.contigs, 'r') as contigs:
        for i, line in enumerate(contigs):
            if line[1:-1] in human_contigs:
                line = contigs.readline()
            else:
                non_human_contigs += line
                
    with open(args.clean, 'w') as clean_contigs:
        clean_contigs.write(non_human_contigs)