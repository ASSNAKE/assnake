from pkg_resources import iter_entry_points 
from assnake.core.config import read_assnake_instance_config
import os, glob, importlib
from assnake.utils.general import read_yaml

class SnakeModule:
    name = ''
    install_dir = ''
    snakefiles = []
    invocation_commands = []
    initialization_commands = []
    wc_configs = []
    initialization_commands = []
    results = []
    dataset_methods = {}
    vizualisation_methods = {}

    def __init__(self, name, install_dir, snakefiles = [], invocation_commands = [], initialization_commands = [], wc_configs = [], results = [], dataset_methods = {}):

        self.assnake_config = read_assnake_instance_config()

        self.name = name
        self.install_dir = install_dir

        results_in_module = glob.glob(os.path.join(self.install_dir, '*/result.py'))
        results_in_module = [ '.'.join(m.split('/')[-3:])[0:-3] for m in results_in_module]
        results_in_module = [ getattr(importlib.import_module(m), 'result') for m in results_in_module ]

        # read_default_config()


        self.results = results
        self.results += (results_in_module)

        self.module_config = self.read_deployed_config()
        print(self.module_config)

        self.snakefiles = snakefiles
        self.invocation_commands = invocation_commands
        self.initialization_commands = initialization_commands
        self.wc_configs = wc_configs
        self.dataset_methods = dataset_methods

    def deploy_module(self):
        if self.assnake_config is not None:
            for result in self.results:
                if result.preset_manager is not None:
                    result.preset_manager.deploy_into_database()

    def read_deployed_config(self):
        if self.assnake_config is not None:
            def_loc = os.path.join(self.assnake_config['assnake_db'], 'module_configs', self.name + '.yaml')
            if os.path.isfile(def_loc):
                return read_yaml(def_loc)
        return None


    @staticmethod
    def get_all_modules_as_dict():
        # Discover plugins
        discovered_plugins = {
            entry_point.name: entry_point.load()
            for entry_point in iter_entry_points('assnake.plugins')
        }

        return discovered_plugins

    