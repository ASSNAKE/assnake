import os
import glob

location = '/data5/bio/databases/fna/handplaced/RNAbacteriophages/*/*.fa'

ff = glob.glob(location)
with open ('./res.fa', 'w') as res:
    for f in ff:
        with open(f,'r') as fasta:
            for line in fasta:
                if line[0] == '>':
                    header = '>' + f.split('/')[-2] + '__' + line[1:-1] + '\n'
                    print(header)
                    res.write(header)
                else:
                    res.write(line)

import subprocess

completed = subprocess.run("""awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' res.fa > result.fa""", shell=True)
print('returncode:', completed.returncode)

