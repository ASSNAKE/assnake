import glob
import os
import sqlite3
import datetime
import yaml
import pkg_resources


wc_config_loc = os.path.join(config['assnake_install_dir'], 'wc_config.yaml')
wc_config = {}
with open(wc_config_loc, 'r') as stream:
    try:
        wc_config = yaml.load(stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print(exc)

discovered_plugins = {
    entry_point.name: entry_point.load()
    for entry_point in pkg_resources.iter_entry_points('assnake.plugins')
}

for module_name, module_class in discovered_plugins.items():
    include: os.path.join(module_class.install_dir, module_class.snakefiles[0])

