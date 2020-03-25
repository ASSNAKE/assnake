import glob, os, pkg_resources, time
import assnake.utils 

wc_config = assnake.utils.load_wc_config()

start = time.time()

# Discover plugins
discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point in pkg_resources.iter_entry_points('assnake.plugins')
}

# We need to update wc_config first
for module_name, module_class in discovered_plugins.items():
    config.update({module_name:module_class.install_dir})
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
        
    for res in module_class.results:
        for sn in res.workflows:
            include: os.path.normpath(sn)
            pass
