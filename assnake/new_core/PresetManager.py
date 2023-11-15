import click
from typing import List

import os, shutil, json, glob, yaml
import zlib

from assnake.new_core.Dataset import Dataset

class Preset:
    def __init__(self, file_path, dataset_name=None):
        self.file_path = file_path
        self.full_name = os.path.basename(file_path)
        
        self.dataset_name = dataset_name
        self.is_global = dataset_name is None
       
        # Splitting file name to extract name, hash, and format
        parts = self.full_name.split('.')
        self.name = parts[0]
        self.current_hash = parts[1] if len(parts) == 3 else None
        self.format = parts[-1]

        self.name_wo_ext = f"{self.name}.{self.current_hash}"
        # Check if the preset is updated and perform update if needed
        if not self.is_updated():
            self.update_filename()

        
    def compute_hash(self):
        with open(self.file_path, 'rb') as f:
            file_content = f.read()
            hash_crc32 = zlib.crc32(file_content)
            return format(hash_crc32, 'x')

    def update_filename(self):
        new_hash = self.compute_hash()
        dir_name = os.path.dirname(self.file_path)
        new_name = f"{self.name}.{new_hash}.{self.format}"
        new_file_path = os.path.join(dir_name, new_name)
        
        if new_file_path != self.file_path:
            os.rename(self.file_path, new_file_path)
            self.file_path = new_file_path
            self.full_name = new_name
            self.current_hash = new_hash
            self.name_wo_ext = f"{self.name}.{self.current_hash}"

    def copy_to_location(self, destination_dir):
        destination_path = os.path.join(destination_dir, os.path.basename(self.file_path))
        os.makedirs(destination_dir, exist_ok=True)
        shutil.copyfile(self.file_path, destination_path)
        return destination_path

    def get_contents(self):
        with open(self.file_path, 'r') as file:
            if self.format == 'json':
                return json.load(file)
            elif self.format == 'yaml':
                return yaml.safe_load(file)
            elif self.format == 'txt':
                return file.read()
            else:
                raise ValueError("Unsupported file format.")

    def is_updated(self):
        return self.compute_hash() == self.current_hash

    def to_cli_option(self):
        def decorator(func):
            return click.option('--' + self.name, default=self.file_path, help=f'Preset for {self.name}', type=click.Path())(func)
        return decorator


class PresetManager:
    def __init__(self, dir_in_database, included_presets_dir, static_files_dir_name, preset_file_format, module_name, result_name):
        self.dir_in_database = dir_in_database
        self.included_presets_dir = included_presets_dir
        self.static_files_dir_name = static_files_dir_name
        self.preset_file_format = preset_file_format
        self.module_name = module_name
        self.result_name = result_name
        self.presets = []

    def load_presets(self, dataset_name=None):
        global_presets = self._find_presets_in_directory(self.dir_in_database)
        dataset_specific_presets = self._find_presets_in_directory(self._dataset_preset_dir(dataset_name)) if dataset_name else []
        self.presets = global_presets + dataset_specific_presets
    
    def find_preset_by_name(self, preset_name, dataset=None):
        self.load_presets(dataset)
        return next((preset for preset in self.presets if preset.name == preset_name and preset.dataset_name == dataset), None)
    
    def find_preset_by_name_and_dataset(self, preset_name, dataset_name):
        """
        Searches for a preset by name within a specific dataset. If not found in the dataset
        but exists in the global directory, copies it to the dataset directory.

        Args:
            preset_name (str): The name of the preset to search for.
            dataset_name (str): The name of the dataset where the preset is searched.

        Returns:
            Preset: The found or copied preset, or None if not found.
        """
        # Load presets from the dataset directory
        self.load_presets(dataset_name)

        # Try to find the preset in the loaded presets
        preset = next((p for p in self.presets if p.name == preset_name and p.dataset_name == dataset_name), None)

        # If not found in the dataset, try the global directory
        if preset is None:
            self.load_presets()  # Load global presets
            global_preset = next((p for p in self.presets if p.name == preset_name and p.is_global), None)

            if global_preset:
                # Copy the global preset to the dataset directory
                copied_preset_path = self.copy_preset_to_dataset(global_preset, dataset_name)
                copied_preset = Preset(copied_preset_path, dataset_name=dataset_name)
                
                # Add the copied preset to the collection and return it
                self.presets.append(copied_preset)
                return copied_preset

        return preset  # Return the found preset or None

    def copy_preset_to_dataset(self, preset, dataset_name):
        dataset = Dataset(dataset_name)
        if preset:
            return preset.copy_to_location(os.path.join(dataset.full_path, 'presets', self.result_name))

    def deploy_into_database(self):
        for preset_file in glob.glob(os.path.join(self.included_presets_dir, f'*.{self.preset_file_format}')):
            preset = Preset(preset_file)
            if not preset.is_updated():
                preset.update_filename()
            preset.copy_to_location(os.path.join(self.dir_in_database, self.result_name))

    def gen_click_option(self):
        self.load_presets()

        if len(self.presets) > 0:
            help_m = f"Preset to use. Available global presets: {[preset.name for preset in self.presets]}"
            help_m += "\n You can view Dataset specific preset by specifying --dataset option and calling help."
            default = self.presets[0]
        else:
            help_m  = 'No presets in database!'
            default = 'No presets in database!'

        return [click.option('--preset',
                        help=help_m,
                        required=False,
                        default = default)]
    
    def _find_presets_in_directory(self, dir_path):
        preset_files = glob.glob(os.path.join(dir_path, f'*.{self.preset_file_format}'))
        return [Preset(preset_file) for preset_file in preset_files]
    
    def _dataset_preset_dir(self, dataset_name):
        dataset = Dataset(dataset_name)
        return os.path.join(dataset.full_path, 'presets') if dataset else ''


