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

## How it works

Assnake assumes that you work with *datasets*, and inside each *dataset* you have multiple *samples*.

The pipeline itself is built using Snakemake. You define rules, which take files as input an produce output files, executing some code. This easy yet powerful technique allows you to build complex data processing pipelines. Snakemake understands what to do by looking at input/output file paths. 

Consider the following: you have your samples, and you want to assemble them and map on resulting contigs. Your rules will look like this (note that we ask for contigs in mapping rule, and input file for mapping rule is the output file for assembly rule, that's the magic):
```
rule assemble:
    input: 
        r1 = '{sample}_R1.fastq.gz',
        r2 = '{sample}_R2.fastq.gz'
    output:
        contigs = '{sample}_assembly.fasta'
    run:
        shell('megahit -1 {input.r1} -2 {input.r2} -o {output.contigs}')

rule map_on_contigs:
    input: 
        r1 = '{sample}_R1.fastq.gz',
        r2 = '{sample}_R2.fastq.gz',
        contigs = {assembled_sample}_assembly.fasta
    output:
        sam = '{sample}_vs_{assembled_sample}_assembly.sam'
    run:
        shell('bwa {input.r1} {input.r2} > {output.sam')

samples = ['A', 'B', 'C']
rule assemble_and_map:
    input: 'A_vs_A_assembly.sam'
```

What we have in curly braces is called *wildcards*. *wildcards* take some user defined value, and that's how Snakemake knows exactly what to do. In the example above rule `map_on_contigs` we have wildcards `sample` and `assembled_sample`. We can give `sample` value `A` and `assembled_sample` value `B`, thus telling snakemake to assemble `B` and map `A` on resulting contg. *Wildcards* are the central part of Assnake, and understanding them is vital.

## Data structure

Metagenomic data is stored in `{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz`. Let's break down the wildcards:

* `prefix` - denotes absolute directory where your dataset is stored. Usage of prefix allows you to store data on different disks, or structure your data by host species, and generally provides more flexibility.
* `df` - name of your dataset.
* `preproc` - is short for *preprocessing* and carries information about how you preprocessed your data. You can cascade preprocessing steps using `__` as a delimeter. For example, `raw__tmtic_def1__no_hum` means that you processed yor raw data with trimmomatic usind default1 parameters and than removed human DNA contamination.
* `sample` - name of your sample.
* `strand` - strand of your sequences. Assnake works only with paired-end sequences at the moment.

Meta-information is stored inside `assnake_db` directory of your choice. The structure of this folder is as follows:

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

* `df_info.yaml` holds meta-info about dataset. Mandatory fields are:
    * df - name of the dataset. Must be the same as dataset dir name.
    * fs_prefix - absolute path of the directory, where dataset is stored.
* `sources.tsv`
* `biospecimens.tsv`
* `mg_samples.tsv`



### Reference data
Reference data is stored inside the `database` directory. It includes various databases for tools and fasta files. Fasta files are stored in `fna_db_dir` and inside that folder you can nest folders as deep as you like. For example: `handplaced/RNAbacteriophages/Escherichia_virus_Qbeta/Escherichia_virus_Qbeta_genomic.fa` When asking for output mapped file just replace `/` with `__` in path. 

### What can be done with Assnake

#### Prepocessing

* Trimmomatic
    * Resulting file to request: `{prefix}/{df}/reads/{preproc}__tmtic_{params}/{sample}/{sample}_R1.fastq.gz`
* FastQC
    * Resulting file to request: `{prefix}/{df}/reads/{preproc}/{sample}/profile/{sample}_{strand}_fastqc.zip`
* BBmap repair
    * Resulting file to request: `{prefix}/{df}/reads/{preproc}__repair/{sample}/{sample}/_R1.fastq.gz`

#### Taxonomic profiling


### Assembly
You need to configure `assembly_dir` where assembles will be stored. Final contigs will be stored in `fna_db_dir`.

File to request: `{fna_db_dir}/assembly/mh__{params}/{dfs}/{samples}/{preprocs}/final_contigs__{min_len}.fa`

dfs are separated by `+`, preprocs by `--` and samples by `:`.

## Implementing new processing step.
In order to implement new processing step you need to create directory `./results/<result_name>`. Take any of the files Asshole can generate right now, or write your own rules.

Than you have a couple of options:

1. The easy way. Just create a file `./results/<result_name>/<result_name>.py` and code your snakemake rule there. Than include this file into `./bin/snake/base.py`. 
