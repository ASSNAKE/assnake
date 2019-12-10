import os
import pandas as pd
import argparse


# contigs_centr_class
# contigs_fasta_loc
# clean_contigs_fasta_loc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--classification', type=str, default='',
                        help='Centrifuge classification')
    parser.add_argument('--out', type=str, default='',
                        help='File to write list of human contigs')
    args = parser.parse_args()

    ccc = pd.read_csv(args.classification, sep = '\t')

    human_contigs = list(ccc.loc[ccc['taxID'] == 9606]['readID'])

    print('Percent of human contigs: ', len(human_contigs)/len(ccc))
    skip = False

    with open(args.out, 'w') as human_contigs_file:
        for item in human_contigs:
            human_contigs_file.write("{}\n".format(item))
    
        

if __name__ == '__main__':
    main()