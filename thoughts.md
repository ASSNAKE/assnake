So, we have two distinct parts of the system. 
The first one is the assnake itself which works with database (in form of files on the disk to streamline things)
The second one is the snakemake pipeline. Pipeline uses some functions from assnake API in order to build the rules correctly when 
the output files cannot be determined directly from the wildcards given in the input. The most obvious example is the task of assembly.
Files for the assembly are taken from the sample_set.tsv file, where samples from different datasets, but the fs_prefix is not given.

We have a modular structure for the pipeline for the easier addition of new analysis tools. We also have a central command for requesting results. The question is: what is the best strategy for coding request behavior? Id the rule takes one sample as input, no parameters and generates the output everything is straightforward. Get the wildcard string from the json included with the module, put the right values and you are good to go. If the rule need some parameters we can parse the included params schema and ask the user to input all the values, than request the name for this params, put it to the db and we are fine again. What if the user
Well the actual question is: do we need to write custom click commands for every result, or can we make it abstract? Or should we mix these to approaches?
If we choose custom commands:
1. + We can control everything precisely and we can add as many native click arguments to the command as we want, including all the parameters and anything rule specific.
2. + It makes the architecture even more modular 
3. - We loose points in terms of uniformity of the system
4. - How the fuck will we request multiple results in one run? Even something basic as metaphlan2 + fastqc + count? The workaround would be to first construct the snakemake file with all the code we need to get the results and than run it. Sounds like a viable option and this way we can add as many datasets samples and results as we want, no limitations what so ever. Sound like a viable option. This snakemake file can be even modified by hand. But where to store this file? Inside the assnake database dir? In the run folder even maybe? Or can we run multiple commands (chain commands)? Well we can chain commands easily indeed. Shit, the problem is that snakemake is called at the end of the command...
This issue is solved by putting all the results inside the global config variable. At the end of the chained commands we should call that special command that runs everything. How to check the order of the chained commands?

well okay fastqc implemented and it works as expected (pretty much)
The new problem is: where should constructing the file path for result go? In module? In sample set?

Now I will try to create an isolated command for fastqc, load it at runtime live a happy life.
https://stackoverflow.com/questions/54073767/command-line-interface-with-multiple-commands-using-click-add-unspecified-optio
So you run something like:
`assnake result request fastqc -d DATASET -p PREPROCESSING ` 

Also we need to be able to conveniently request results and run snakemake from the web and celery. 

Now the process looks like this:
1. User requests results
2. Output result file paths are constructed
3. Snakemake gets called

Study is the most top level object and can include unlimited datasets of different nature (WGS in 16s rRNA are supported at the moment).
Dataset object
SampleSet object - this object can include arbitrary amount of samples from any datasets, but most commonly ot is used to combine representative samples for one dataset. 
SequencingRun object


# 17.01.2020
Back again to the organisation of the plugin architecture
