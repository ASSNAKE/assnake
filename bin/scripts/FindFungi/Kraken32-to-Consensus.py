"""
Script that takes the concatenated sorted output of the Kraken32
databases and finds a consensus prediction for each read
"""

__author__ = "Paul Donovan"
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovandonegal@gmail.com"

from itertools import islice
import sys
from collections import Counter
import argparse

# Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Sorted.Kraken.tsv')
parser.add_argument('Consensus.Kraken.tsv')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

# Define lists and variables
n = 32  # The number of lines to loop over at once
GuessList = list()
YesCount = 0
NoCount = 0

def find_consensus(next_n_lines, GuessList, ):
    # YesCount = YesCount + 1
    count = Counter(GuessList)  # #
    freq_list = count.values()  # #
    max_cnt = max(freq_list)  # Find the most common item in the list #
    total = freq_list.count(max_cnt)  # #
    most_common = count.most_common(total)  # #
    if len(most_common) == 1:  # If there is one most common taxid
        read = (next_n_lines[0].strip().split("\t")[1])
        return("C\t" + read + "\t" + (most_common[0])[0] + "\t100\t100\n")
    else:  # If there is more than one most common prediction
        prediction_list = list()
        for i in next_n_lines:  # Loop over the 32 kraken predictions again
            item = i.strip().split("\t")
            if str(item[2]) == str(0):  # If kraken failed to predict a taxid for the read (aka 0), skip this
                pass
            else:
                Prediction = item[4].split(" ")  # Gather all predictions from all kraken runs
                for u in Prediction:
                    if "A" not in str(u):  # Don't use the 'ambiguous' flag in our taxonomy predictions
                        prediction_list.append(u)  # Add all predictions to a list
        new_dict = dict()
        for kmer in prediction_list:  #
            taxid = (kmer.split(":"))[0]  #
            kmer_number = (kmer.split(":")[1])  #
            if str(taxid) == str(0):  #
                pass  # Loop over all predictions from all 32 kraken runs
            else:  # and create a dictionary where taxid is the key
                if str(taxid) in new_dict:  # and the number of kmers supporting the taxid is
                    new_dict[taxid] = (int(kmer_number) + int(new_dict[taxid]))  # the value
                else:  #
                    new_dict[taxid] = int(kmer_number)  #
        if bool(new_dict) == False:
            print('this is the case')# If the dictionary is empty, don't do anything. This should never be the case
            pass
        else:
            TaxidPrediction = str(max(new_dict,
                                      key=new_dict.get))  # Get the key with the highest value. Draws may occur but a key will be selected
            TypicalRead = next_n_lines[0].strip().split("\t")
            return ("C\t" + TypicalRead[1] + "\t" + TaxidPrediction + "\t100\t100\n")
        new_dict.clear()
# Open output file for writing
Output = open(sys.argv[2], "w")


# Find the consensus prediction using 32 predictions for the same read
with open(sys.argv[1]) as f:

    current = f.readline()
    while current:
        next_line = f.readline()
        if next_line:
            current_read = current.split("\t")[1]
            next_read = next_line.strip().split("\t")[1]
            GuessList = [(current.strip().split("\t"))[2]]
            next_n_lines = [current]
            while current_read == next_read:
                GuessList.append((next_line.strip().split("\t"))[2])
                next_n_lines.append(next_line)
                next_line = f.readline()
                if next_line:
                    next_read = next_line.strip().split("\t")[1]
                else:
                    next_read = ''
            current = next_line
            Output.write(find_consensus(next_n_lines, GuessList))
        else:
            break

    print "Done"

Output.close()
