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
    parser.add_argument('--contigs', type=str, default='',
                        help='Fasta with contigs')
    parser.add_argument('--clean', type=str, default='',
                        help='File to write clean contigs')
    args = parser.parse_args()

    print(args.classification, args.contigs, args.clean)
    ccc = pd.read_csv(args.classification, sep = '\t')

    human_contigs = list(ccc.loc[ccc['taxID'] == 9606]['readID'])

    print('Percent of human contigs: ', len(human_contigs)/len(ccc))
    non_human_contigs = ''
    skip = False

    with open(args.clean, 'w') as clean_contigs:
        with open(args.contigs, 'r') as contigs:
            for i, line in enumerate(contigs):
                if skip:
                    skip = False 
                else:
                    if line[1:-1] in human_contigs:
                        skip = True
                        #line = contigs.readline()
                    else:
                        clean_contigs.write(line)
                        # non_human_contigs += line
                
    
        

if __name__ == '__main__':
    main()