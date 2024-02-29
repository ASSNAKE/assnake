## Create dataset and add test data

After installing the program, navigate to the directory your created for your Dataset, for example `cd /mnt/disk1/RUNS/your_username/MY_DATASET`

## Register this folder as a dataset in the system 

Inside the folder, run `assnake dataset init` to register the folder as the dataset in the metasnake system. 

## Add test data

This guide is designed to start from the constructed feature table. Obtain the files from the developer on Telegram. Then inside the Dataset folder create new directory that will hold the feature table. For example, `mkdir -p feature_tables/all_samples/dada2_noPool`
In wildcard terms that would be `feature_tables/{sample_set}/{ft_name}`. Put the `phyloseq.rds` and `metadata.yaml` into this folder.

## Install the packages that are designed for downstream analysis 

Navigate to the folder that holds the code (something like `~/MetasnakeGit`), and perform the following:

```bash
git clone https://github.com/ASSNAKE/assnake-phyloseq.git
cd assnake-phyloseq
pip install -e ./
```

You also need to have metasbmR package installed. 
You can obtain it here https://github.com/META-SBM/metasbmR

Clone the repo, then in R install the package, like this 

```R
devtools::document('~/ProjectsGit/metasbmR')
devtools::build('~/ProjectsGit/metasbmR')
devtools::install('~/ProjectsGit/metasbmR')
```


## Create recipe.py file and test that it works

`recipe.py` file defines what Results we want to produce for a given set of data. 


```python
from assnake.core.Pipeline import Pipeline

# Create pipeline instance
pipeline = Pipeline(name = 'testing_tests', description = 'testing')

# Specify what you want
heatmap_params = {
    'dataset': 'MY_DATASET', 
    'feature_table_name': 'dada2_noPool', 
    'sample_set': 'all_samples', 
    # Note that we are not explicitly requesting filtering
    # All the relevant info is passed by this string, metasnake will take care of underlying parsing. 
    'filter_chain': 'filter_taxa_by_detection_1_prevalence_0.03/transform_clr', 

    'heatmap_preset': 'None'
}
# Add your request to pipeline
pipeline.add_result('pheatmap', **heatmap_params)
```


Save this file, preferably in the root of the dataset directory. `/mnt/disk1/RUNS/your_username/MY_DATASET/recipe.py`
Now you can test if pipeline works by running

`assnake pipeline pipeline  --pipeline-script /mnt/disk1/RUNS/your_username/MY_DATASET/recipe.py --run`

If it works - great, if not - contact the developer on Telegram. 