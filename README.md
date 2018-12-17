# assnake
Snakemake pipeline for metagenomics data. Part of Boimedical Research Platform.

# WORK IN PROGRESS. NOT READY FOR PRODUCTION.

## Installation
1. You need a conda environment with Snakemake 5+ installed.
2. Download pipeline and place it into desired directory.
3. Edit config.yml to provide paths to your local tools. 

## Data storage
All paths are relative to pipeline folder.
Reads live in `datasets/{df}/reads/{preproc}/{sample}/`
`preproc` is short for preprocessing. Import your reads to `imp` folder. 
You can preprocess your reads by creating chain preprocessing names, for example `imp__tmtic_def1` means that imported reads were processed with Trimmomatic with def1 parameters.


