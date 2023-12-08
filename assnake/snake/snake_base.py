import glob, os, time
from assnake.core.config import load_wc_config
from pkg_resources import iter_entry_points 
wc_config = load_wc_config()

start = time.time()


def get_previous_step_output(wildcards):
    # Extract the current step number and reconstruct the previous step's output path
    # current_step_num = int(wildcards.step_num)
    print(wildcards)
    input_target = ''
    # if current_step_num > 1:
    #     input_target =  f"{wildcards.fs_prefix}/{wildcards.df}/feature_tables/{wildcards.sample_set}/{wildcards.ft_name}/{wildcards.filter_chain.strip('/')}/phyloseq.rds"
    # else:
    #     # Return the initial phyloseq file path for the first step
    #     input_target =  f"{wildcards.fs_prefix}/{wildcards.df}/feature_tables/{wildcards.sample_set}/{wildcards.ft_name}/phyloseq.rds"

    input_target =  f"{wildcards.fs_prefix}/{wildcards.df}/feature_tables/{wildcards.sample_set}/{wildcards.ft_name}/{wildcards.filter_chain.strip('/')}/phyloseq.rds"

    return(input_target.replace('//', '/'))


# Discover plugins
discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point in iter_entry_points('assnake.plugins')
}


# We need to update wc_config first
for module_name, module_class in discovered_plugins.items():
    
    module_config = {'install_dir': module_class.install_dir}

    if module_class.module_config is not None:
        module_config.update(module_class.module_config)

    config.update({module_name:module_config})


    for wc_conf in module_class.wc_configs:
        if wc_conf is not None:
            wc_config.update(wc_conf)
    for res in module_class.results:
        if res.wc_config is not None:
            wc_config.update(res.wc_config)


# and now include all the stuff
for module_name, module_class in discovered_plugins.items():

    for snakefile in module_class.snakefiles:
        include: os.path.normpath(os.path.join(module_class.install_dir, snakefile))
        print(os.path.normpath(os.path.join(module_class.install_dir, snakefile)))

    for res in module_class.results:
        # print(res)

        for sn in res.workflows:
            include: os.path.normpath(sn)
