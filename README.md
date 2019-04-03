# ASSNAKE

Assnake is a pipeline for metagenomics data analysis. It is built with Snakemake and is fully modular. It tracks all your tool versions, parameters and environments.


Analysis of metagenomics data consists of 2 steps: quality control and preprocessing, and analysis itself. 

## List of implemented preprocessing steps:
1. Downloading data from SRA.  
2. Count number of base pairs and reads in fastq files.
3. Check quality with FastQC.
4. Trimming with Trimmomatic.

## List of implemented analysis steps and tools:
* Taxonomic profiling.
    1. Metaphlan2
    2. Kraken 1
    3. Centrifuge
* Functional profiling
    1. Humann2
* Mapping
    1. Bowtie2
    2. BWA
* Assembly
    1. Megahit

## Architecture.
Assnake assumes that you work with *datasets*, and inside each *dataset* you have multiple *samples*. All *dataset* names inside of assnake instance **must** be unique. All *sample* names inside one *dataset* **must** be unique.

Meta-information is stored inside `database` directory of your choice. The structure of this folder is as follows:

```
|-- <df_name_1>
   |-- df_info.yaml
   |-- sources.tsv
   |-- biospecimens.tsv
   |-- mg_samples.tsv
|-- <df_name_2>
   |-- df_info.yaml
   |-- sources.tsv
   |-- biospecimens.tsv
   |-- mg_samples.tsv
...
```

`df_info.yaml` holds meta-info about dataset. Mandatory fields are:
1. df - name of the dataset. Must be the same as dataset dir name.
2. fs_prefix - prefix of the dataset location on the file system. 

Metagenomics data is stored in `{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz`. Usage of prefix allows you to store data on different disks, or structure your data by host species, and generally provides more flexibility. `preproc` is short for *preprocessing* and carries information about how you, well, *preprocessed* your data :). You can cascade preprocessing steps using `__` as a delimeter. 

### Reference data
Reference data is stored inside the `database` directory. It includes various databases and fasta files. Fasta files are stored in `fna_db_dir` and inside that folder you can nest folders as deep as you like. For example: `handplaced/RNAbacteriophages/Escherichia_virus_Qbeta/Escherichia_virus_Qbeta_genomic.fa`

### Assembly
We allow assembly from any combinations of samples. You need to configure `assembly_dir`where assembles will be stored.
Final contigs are stored in `fna_db_dir`.

## Implementing new processing step.
In order to implement new processing step you need to create directory `./results/<result_name>`.
Than you have a couple of options:
1. The easy way. Just create a file `./results/<result_name>/<result_name>.py` and code your snakemake rule there. Than include this file into `./bin/snake/base.py`