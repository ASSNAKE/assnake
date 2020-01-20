from snakemake.shell import shell


shell('''rm -rf {snakemake.params.out_folder};\
     spades.py --meta -1 {snakemake.input.r1} -2 {snakemake.input.r2} -o {snakemake.params.out_folder} -t {snakemake.threads} > {snakemake.log} 2>&1;\
         cp {snakemake.params.contigs} {snakemake.output.final_contigs}''')