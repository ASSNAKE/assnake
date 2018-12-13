# !/home/fedorov/miniconda3/bin/python

# Simple wrapper to trim reads in gzipped fastq file to desired length and save to new fastq.gz file
import gzip
from Bio import SeqIO
import argparse
import os.path

def read_args():
    # Read input arguments
    parser = argparse.ArgumentParser()
    
    g1 = parser.add_argument_group('Input options')
    g1.add_argument('--input', help='Input FASTQ file', required=True)
    g1.add_argument('--out', help='Output FASTQ file', required=True)
    g1.add_argument('--len', help='Leave only `len` reads', type=int, default=100)
    
    return parser.parse_args()

def trim_reads_in_file(in_fastq_loc, out_fastq_loc, len_to_trim = 100):
    handle_out = gzip.open(out_fastq_loc, "wt")
    print('Opened ' + out_fastq_loc + ' for writing')
    
    with gzip.open(in_fastq_loc, "rt") as handle:
        i = 0
        print('Opened ' + in_fastq_loc + ' for reading')
        for record in SeqIO.parse(handle, "fastq"):
            if len(record.seq) > len_to_trim:
                handle_out.write(record[0:len_to_trim].format("fastq"))
            else:
                handle_out.write(record.format("fastq"))
            if i % 100000 == 0:
                print('INFO: Finished with read number: ' + str(i))
            i += 1
    handle_out.close()
            
def main():
    args = read_args()
    if os.path.isfile(args.input) and args.len > 1:
        print('Starting trimming')
        trim_reads_in_file(args.input, args.out, args.len)
    else:
        print("ERROR")
    
if __name__ == "__main__":
    main()