__author__ = "Dmitry Fedorov"
__copyright__ = "Copyright 2019, Dmitry Fedorov"
__email__ = "fedorov.de@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


shell('''metaphlan2.py --mpa_pkl {snakemake.params.MPA_PKL} --bowtie2db {snakemake.params.BOWTIE2DB} \
         {snakemake.input.r1},{snakemake.input.r2} --input_type fastq --bowtie2out {snakemake.params.b} \
         --nproc {snakemake.threads} > {snakemake.output.o}''' )