# ASSNAKE

Assnake is a system for NGS data analysis, vizualisation and management.
It allows you to go from raw reads to biological insights in just a few commands.

As of pre-alpha release Assnake is capable of full-blown metagenomic data analysis (Both Shotgun WGS and Amplicon 16s included!) and some RNA-seq analysis.

Assnake was born in an effort to provide userfriendly, scalable and powerful system for NGS data analysis, but at the same time our goal was to make it easy to extend with your own pipelines.

Assnake has 3 key concepts:

* You store reads of your Samples inside Datasets.
    Datasets are just folders on your file system, assnake assumes that raw reads are stored inside `{PATH_TO_DATASET_FOLDER}/{YOUR_DATASET_NAME}/reads/raw`
    So, for example you may have your reads stored at /home/ozzy/bat_microbiome/reads/raw. /home/ozzy is the prefix of your Dataset and bat_microbiome is it's name. You need to register your dataset inside Assnake with `assnake dataset create -f /home/ozzy -d bat_microbiome` Don't be afraid, Assnake is very careful and will never overwrite or delete your data!
* Data needs quality control and cleaning before being analyzed.
* Everything that produces meaningful and useful data produces Result which you can request from Assnake `assnake result <RESULT_NAME> --df <DATASET> run`. For example command `assnake result fastqc --df bat_microbiome run` will produce fastqc reports for all samples in Dataset bat_microbiome.


Assnake draws inspiration mainly from Anvio (incorporated into Assnake) and QIIME-2. 

Assnake uses Snakemake as a workflow subsystem, thus it gets Snakemake's ability to run on all kinds of servers, clusters, personal computers and in the cloud.

Assnake is entirely open sourced and is built with great open-source community. We encourage the community to join Assnake's initiative in creating open reproducable and userfriendly omics data analysis and management.

# Quick start
## Conda part
1. Install conda https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
2. Create new environment `conda create -n assnake python=3.6`
3. Activate your new environment `source activate assnake`
## Assnake part
1. Navigate to some directory on your file system, for example your home directory `cd ~`
2. Clone this repository to your computer using `git clone https://github.com/Fedorov113/assnake.git` 
3. Enter to just created assnake directory `cd assnake`
4. Install the package using `pip install -e ./`
5. Verify your installation by running `assnake --help`


## Initialization
Call `assnake init start`

It will ask you which directory you would like to use for assnake database, choose some folder on file system with at least 20 Gb of free space. If the folder is not yet created, assnake will create it. 

Now you can start processing your data!

Now we need to register or create dataset in assnake. Run `assnake dataset create --df <DATASET_NAME> --fs_prefix <FOLDER WHERE TO STORE DATA>`
Full folder name looks like this: '{fs_prefix}/{df}'. Say, I want to create dataset with name `miseq_sop` and I want to put in `/home/fedorov/bio`.
I will call `assnake dataset create --df miseq_sop --fs_prefix /home/fedorov/bio`. Folder `/home/fedorov/bio/miseq_sop` will be created on file system, if not created already. If you already have same folder it is totally OK.

Now run `assnake dataset info -d <DATASET_NAME>`

## Installation of modules
This version is for Vlad and dada2, so we need only 
1. https://github.com/Fedorov113/assnake-dada2
2. https://github.com/Fedorov113/assnake-core-preprocessing 

Just clone this repositories and install with `pip install -e ./` while `assnake` conda env is activated (`source activate assnake`)

After installation run `assnake result request --help` and you will see new available results.

# Running DADA2 for 16s rRNA data
Run command without `< >` symbols around your dataset name.
1. Run `assnake result request dada2-filter-and-trim -d <YOUR_DATASET> -p raw run --threads 1 --jobs 4 --run`. This will filter your reads by quality with default parameters using 4 jobs in parallel and 1 thread on each job.
2. Execute `assnake result request dada2-full -d <YOUR_DATASET> -p raw__dada2fat_def run -t 4 -j 1 --run`. Now we run 1 jobs with 4 threads. If you know that your machine has more available cores, feel free to use them and increase threads or jobs. 

You are done! You can find dada2 results at `{FS_PREFIX}/{YOUR_DF}/dada2/sample_set/learn_erros__def/seqtab_nochim__20.rds` and `{FS_PREFIX}/{YOUR_DF}/dada2/sample_set/learn_erros__def/taxa_20.rds`.
Just load this files in R using `readRDS()` function.


<!-- # OLD

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
* Antibiotic resistance annotation.
    1. Megares

## Configuration

You need to configure the following parameters in your `config.yml`

```
assembly_dir: 'absolute path to direcrory for storing assemblies'
assnake_db: 'absolute path to directory for storing meta-information'
fna_db_dir: 'absolute path to your fasta database'
```

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
* `sources.tsv` contains information about objects of study (sources of biospecimens). Columns are:
    * `source` - unique identificator of source. Alphanumeric, `-` and `_`
    * `description` - free text with description of source.
* `biospecimens.tsv`
    * `biospecimen`
    * `source`
    * `time`
    * `description`
* `mg_samples.tsv`
    * `biospecimen`
    * `fs_name`



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

## Metagenome Assembled Genomes
MAG workflow is built around anvio, which is an awesome tool with great capabilities, you should defenitley check it out.
MAGs are binned using contigs assembled fron *n* samples. First, we create contigs_db in the folder with contigs, then we map reads back to contigs, create profiles and after that we apply binning step. #TODO We should store sample_profiles in samples mapping folder. We store merged profile in `{prefix}/{df}/anvio/merged_profiles/bwa__{bwa_params}___assembly___mh__{params}___{dfs}___{samples}___{preprocs}/MERGED/`



## Implementing new processing step.
In order to implement new processing step you need to create directory `./results/<result_name>`. Take any of the files Asshole can generate right now, or write your own rules.

Than you have a couple of options:

1. The easy way. Just create a file `./results/<result_name>/<result_name>.py` and code your snakemake rule there. Than include this file into `./bin/snake/base.py`. 

# Running tests
`python cli.py -s ./tests/test_snake.py`
 -->
