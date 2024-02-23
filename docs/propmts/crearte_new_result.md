# How to create new Result?

## In order to create new Result in Assnake you have to do the following:

1. Decide what program or script you want to wrap as an Assnake result. If you have a chunk of code in, say R, create a function for it, add parameters you want to control as arguments. 
2. Prepare wrapper method for a function that accepts input and output filenames, reads the input calls computing function and writes the result to destination. (optional)
3. Write a snakemakes workflow.smk file defining the rules for the result, call your script or program in shell directive.
4. Create a result.py file, where you create assnake Resukt object, specifying name, description, input_type and other attributes.
5. Create wc_config.yaml file, where you specify `{result_name}_wc` that will be used by Assnake for target file generation. 


## Here is a concrete example of this process:

### Extracting the script

We have this chunk of code for plotting bar charts:

```R
p <- ps_all_cov %>%
  comp_barplot(
    tax_level = "Genus", 
    n_taxa = 20, 
    other_name = "Other", 
    # tax_transform_for_ordering = 'compositional', 
    order_with_all_taxa = TRUE,
    sample_order = "bray",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 20, add = "grey90"),
    merge_other = FALSE, 
    bar_outline_colour = NA, # Set to NA to remove grey boxes
    bar_width = 1
  ) +
  labs(y = NULL, x =NULL) +
  facet_wrap("batch_x", nrow = 4, scales = "free") +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 75, hjust = 1) # Rotate x labels
  )

ggsave("bar_facet_v0.pdf", plot = p, device = "pdf", width = 300, height = 135, limitsize = FALSE)
```

We create a wrapper function from it and save it to:
`./barplot_facet/rscript.R`

```R
plot_bar_facet_v0 <- function(filepath, output_pdf, n_taxa = 20, tax_level = "Genus", facet_variable = "batch_x", pdf_width = 300, pdf_height = 135) {
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(microViz)

  ps_obj <- readRDS(filepath)


  p <- ps_obj %>%
    comp_barplot(
      tax_level = tax_level, 
      n_taxa = n_taxa, 
      other_name = "Other", 
      order_with_all_taxa = TRUE,
      sample_order = "bray",
      taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
      palette = distinct_palette(n = n_taxa, add = "grey90"),
      merge_other = FALSE, 
      bar_outline_colour = NA,
      bar_width = 1
    ) +
    labs(y = NULL, x =NULL) +
    facet_wrap(facet_variable, nrow = 4, scales = "free") +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 75, hjust = 1)
    )

  ggsave(output_pdf, plot = p, device = "pdf", width = pdf_width, height = pdf_height, limitsize = FALSE)
}

```


Then we create workflow.smk file:
`./barplot_facet/workflow.smk`

```python
plot_bar_facet_v0_script = os.path.join(config['assnake-phyloseq']['install_dir'], 'barplot_facet/rscript.R')

rule plot_bar_facet_v0:
    input:
        ps = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/phyloseq.rds'
    output:
        pdf = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/bar_facet_v0_tax{tax_level_preset}_n{n_taxa_preset}_facet{facet_variable_preset}_w{pdf_width_preset}_h{pdf_height_preset}.pdf'
    wildcard_constraints:
        df = "[^/]+",
        ft_name = "[^/]+",
        sample_set = "[^/]+",
        tax_level_preset = "[\w]+",
        n_taxa_preset = "\d+",
        facet_variable_preset = "[\w]+",
        pdf_width_preset = "\d+",
        pdf_height_preset = "\d+"
    params:
        plot_bar_facet_v0_script = plot_bar_facet_v0_script
    shell:
        """
        Rscript -e " \
        source('{params.plot_bar_facet_v0_script}'); \
        plot_bar_facet_v0('{input.ps}', '{output.pdf}', {wildcards.n_taxa_preset}, '{wildcards.tax_level_preset}', '{wildcards.facet_variable_preset}', {wildcards.pdf_width_preset}, {wildcards.pdf_height_preset}); \
        "
        """
```

Now for `./barplot_facet/result.py`

```python
import os
from assnake.core.Result import Result

result = Result.from_location(
    name='barplot_facet', # Note that the name must be the same as the directory name!
    description='Plot a facetted barplot.',
    result_type='plot',
    input_type='phyloseq',
    with_presets=False,
    preset_file_format='yaml',
    location=os.path.dirname(os.path.abspath(__file__)))
```

And finally wc_config.yaml `./barplot_facet/wc_config.yaml`
Note that this wildcard strings are the same as in rule plot_bar_facet_v0 output. 
And the the name of the string is `result_name + '_wc'`
```yaml
barplot_facet_wc: '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/bar_facet_v0_tax{tax_level_preset}_n{n_taxa_preset}_facet{facet_variable_preset}_w{pdf_width_preset}_h{pdf_height_preset}.pdf'
barplot_facet_source_wc: '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/bar_facet_v0_tax{tax_level_preset}_n{n_taxa_preset}_facet{facet_variable_preset}_w{pdf_width_preset}_h{pdf_height_preset}.pdf'
```


Great! Now you are ready for integration and testing!


Sometimes we want to pass parameters inside yaml file. In that case in addition to result.py, workflow.smk and wc_config.yaml we also need to create param_schema.yaml

Here is and example implementation of deseq2_plot

==> param_schema.yaml <==
```yaml
formula:
  type: string
  description: The model formula for DESeq2 analysis.
  default: sex + batch

contrast_condition:
  type: string
  description: The condition for the DESeq2 contrast.
  default: condition

contrast_level1:
  type: string
  description: The first level for the DESeq2 contrast.
  default: treated

contrast_level2:
  type: string
  description: The second level for the DESeq2 contrast.
  default: untreated

alpha:
  type: float
  description: Significance level for determining differentially expressed taxa.
  default: 0.05

minlfc:
  type: float
  description: Minimum log fold change to consider a taxon as differentially expressed.
  default: 1
```

==> result.py <==
```python
import os
from assnake.core.Result import Result

result = Result.from_location(
    name='deseq2_plot',
    description='DESeq2 analysis plot.',
    result_type='analysis',
    input_type='phyloseq',
    with_presets=True,
    preset_file_format='yaml',
    location=os.path.dirname(os.path.abspath(__file__)))
```

==> rscript.R <==
```R
save_deseq2_results <- function(ps_filepath, dds_out_file, params_yaml_path, sigtab_out_file, pdf_out_file) {
    params <- yaml::read_yaml(params_yaml_path)
    params

    library(phyloseq)
    library(DESeq2)
    library(ggplot2)
    library(ggpubr)
    library(metasbm)


    ps_obj <- readRDS(ps_filepath)


    results <- plot_deseq2(
        ps_obj,
        formula=params$formula,
        contrast=c(params$contrast_condition, params$contrast_level1, params$contrast_level2),
        alpha=params$alpha,
        minlfc=params$minlfc
    )

    # Save the results
    saveRDS(results[[1]], dds_out_file)
    write.csv(results[[2]], sigtab_out_file, row.names = FALSE)
    ggsave(pdf_out_file, plot = results[[3]], device = "pdf")
}
```
==> wc_config.yaml <==
```yaml
deseq2_plot_wc: '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/deseq2_{preset}_plot.pdf'
deseq2_plot_source_wc: '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/phyloseq.rds'
```

==> workflow.smk <==
```python
deseq2_script = os.path.join(config['assnake-phyloseq']['install_dir'], 'deseq2_plot/rscript.R')


rule deseq2_analysis:
    input:
        ps = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/phyloseq.rds',
        params_yaml = "{fs_prefix}/{df}/presets/deseq2_plot/{preset}.yaml"

    output:
        diagdds = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/deseq2_{preset}_diagdds.rds',
        sigtab = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/deseq2_{preset}_sigtab.csv',
        plot = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/deseq2_{preset}_plot.pdf',
        preset = '{fs_prefix}/{df}/feature_tables/{sample_set}/{ft_name}/{filter_chain}/deseq2_plot_{preset}.yaml'
    wildcard_constraints:
        df = "[^/]+",
        ft_name = "[^/]+",
        sample_set = "[^/]+",
    params:
        deseq2_script = deseq2_script
    shell:
        """
        cp {input.params_yaml} {output.preset}; \
        Rscript -e " \
        source('{params.deseq2_script}'); \
        save_deseq2_results('{input.ps}', '{output.diagdds}', '{input.params_yaml}', '{output.sigtab}', '{output.plot}'); "\
        """
```