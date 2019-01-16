import glob
import os
import sqlite3
import datetime


db_loc = '/data6/bio/TFM/asshole/db.sqlite3'

def save_to_db(task_id, rule_name, in_list, out_list, status ):
    
    save_str_wc = "INSERT INTO explorer_snakeruleresult VALUES (null, '{date_time}', '{task_id}', '{rule_name}', '{in_list}','{out_list}', '{status}');"
    save_str = save_str_wc.format(date_time=datetime.datetime.now(),
task_id=task_id, rule_name=rule_name, in_list=in_list, out_list=out_list, status=status)
    
    print(save_str)
    
    conn = sqlite3.connect(db_loc)
    c = conn.cursor()
    c.execute(save_str)
    conn.commit()
    conn.close()

configfile: 'config.yml'
    
snakefiles = os.path.join(config["software"]["snakemake_folder"], "bin/snake/")
include: snakefiles + "bowtie2.py"
include: snakefiles + "anvio.py"
include: snakefiles + "prokka.py"
include: snakefiles + "bwa.py"
include: snakefiles + "megahit.py"
include: snakefiles + "general.py"
include: snakefiles + "preprocess.py"
include: snakefiles + "marker"
include: snakefiles + "taxa.py"
include: snakefiles + "strain_finder.py"
include: snakefiles + "download.py"
include: snakefiles + "find_fungi.py"
include: snakefiles + "subread.py"



samples = [s.split('/')[-1] for s in glob('/data6/bio/TFM/pipeline/datasets/HELIC/reads/imp/*')]
samples = [s.split('/')[-1] for s in glob('/data6/bio/TFM/pipeline/datasets/DAVID/reads/imp/*B5')]

rule david_count:
    input: expand('datasets/DAVID/reads/imp/{sample}/profile/{sample}_{strand}.count', sample = samples, strand=['R1', 'R2'])

rule david_fastqc:
    input: expand('datasets/DAVID/reads/imp/{sample}/profile/{sample}_{strand}_profile.tsv', sample = ['1712434_B5'], strand=['R1', 'R2'])
#samples = []
# samples = ['DFM_003_F1_S10', 'SFM-017-F1-1_S57', 'T22T5_L1S1_B6_S105']

what_i_want_to_map_in_helic = [
    {'stamm': '1352',
     'ref' : 'A45'
    },
    {'stamm': '9192',
     'ref' : 'A45'
    },
    {'stamm': 'E-14',
     'ref' : 'A45'
    },
    {'stamm': 'hpy',
     'ref' : 'A45'
    },
    {'stamm': '26695',
     'ref' : '26695'
    },
    {'stamm': 'A45',
     'ref' : 'A45'
    },
    {'stamm': 'E-48',
     'ref' : 'E48'
    },
    {'stamm': 'Hab13',
     'ref' : 'H13-1'
    },
    {'stamm': 'J99',
     'ref' : 'J99'
    },
]

def generate_file_loc(info, fasta, ext, gff = False):
    location = []
    wc_str = 'datasets/HELIC/mapped/imp__tmtic_helic1__bbduk_rmaaa/bwa__def/handplaced__helic__{ref}/{ref}_{fasta}/mapped/{stamm}-{rep}'+ext
    if gff:
        fasta_wc = '/data5/bio/databases/fna/handplaced/helic/{ref}/{ref}_{fasta}.gff'
    else:
        fasta_wc = '/data5/bio/databases/fna/handplaced/helic/{ref}/{ref}_{fasta}.fa'

    for element in info:
        fasta_loc = fasta_wc.format(ref = element['ref'], fasta=fasta)
        if os.path.isfile(fasta_loc):
            for i in range(1,4):
                location.append(wc_str.format(rep=str(i), fasta=fasta, **element))
        else:
            pass
            #print(fasta_loc)
    return location

rule helic_map:
    input: generate_file_loc(what_i_want_to_map_in_helic, 'cds_from_genomic', '.bb_stats')

rule helic_fc:
    input: generate_file_loc(what_i_want_to_map_in_helic, 'genomic', '_feature_counts.txt', True)

rule helic_bbduk:
    input: expand('datasets/HELIC/reads/imp__tmtic_helic1__bbduk_rmaaa/{sample}/{sample}_R1.fastq.gz', sample=samples)

rule all:
    input: expand('datasets/FHM/taxa/imp__tmtic_def1/FindFungi/{sample}/all_classified_sorted.tsv', sample = samples)
rule one:
    input: 'datasets/FHM/taxa/imp__tmtic_def1/FindFungi/DFM_003_F1_S10/lca.csv'
rule all_wer:
    input: 'datasets/Tutorial/taxa/imp/FindFungi/p136C/lca.csv'

rule test:
    input: 'datasets/Fungi/taxa/imp/FindFungi/ERR675624/lca.csv'

rule map_helic:
    input: expand('datasets/HELIC/mapped/imp__tmtic_helic1/bwa__def/handplaced__helic__A45/A45_cds_from_genomic/mapped/{sample}.bb_stats', sample=samples)

rule mp2_test:
    input: expand('datasets/Smillie18/taxa/imp/mp2__def/{sample}/{sample}.mp2', sample=['ERR2198714','ERR2198715', 'ERR2198716', 'ERR2198717', 'ERR2198718'])


samples = [s.split('/')[-1] for s in glob('/data6/bio/TFM/pipeline/datasets/RNAViromeFMS/reads/dwn/*')]

rule rna_vir_dwn:
    input: expand('datasets/RNAViromeFMS/reads/dwn__repair/{sample}/{sample}_R1.fastq.gz', sample = samples)

rule _virome:
    input: expand('datasets/RNAViromeFMS/mapped/dwn__repair/bwa__def/handplaced__RNAbacteriophages__combined/combined_v1/mapped/{sample}.bb_stats', sample = samples)


rule rna_virome:
    input: expand('datasets/RNAViromeFMS/mapped/dwn__repair/bwa__def/handplaced__RNAbacteriophages__Pseudomonas_virus_phi6/Pseudomonas_virus_phi6_genomic/mapped/{sample}.bb_stats', sample = samples)
