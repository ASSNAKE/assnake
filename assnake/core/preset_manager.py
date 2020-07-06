import os, glob, click
from assnake.core.config import read_assnake_instance_config
from assnake.utils.general import compute_crc32_of_dumped_dict
import shutil

class PresetManager:

    params_schema = {}
    dir_in_database = ''

    def __init__(self, dir_in_database, included_presets_dir, preset_file_format='json', static_files_dir_name = 'static'):
        self.dir_in_database = dir_in_database
        self.included_presets_dir = included_presets_dir
        self.preset_file_format = preset_file_format
        self.static_files_dir_name = static_files_dir_name

        instance_config = read_assnake_instance_config()
        if instance_config is not None:
            self.preset_file_in_db_wc = os.path.join(
                instance_config['assnake_db'], self.dir_in_database, '{preset}.{hash}.' + self.preset_file_format)

    def deploy_into_database(self):
        # Check that it assnake initialized
        instance_config = read_assnake_instance_config()
        if instance_config is not None:
            # Check for directory for params in current database and create if not
            os.makedirs(os.path.join(
                instance_config['assnake_db'], self.dir_in_database), exist_ok=True)

            # Now get params files we want to import into database
            presets_included_files = glob.glob(
                os.path.join(self.included_presets_dir, '*.' + self.preset_file_format))

            for preset_file_loc in presets_included_files:
                preset_name = os.path.basename(preset_file_loc).split('.')[0]
                preset_crc32_hex = compute_crc32_of_dumped_dict(
                    preset_file_loc)
                # Try to copy default parameters json
                loc_in_db = os.path.join(self.preset_file_in_db_wc.format(
                    hash=preset_crc32_hex,
                    preset=preset_name))
                shutil.copyfile(preset_file_loc, loc_in_db)

            # Now go with the static files
            static_files_locs = glob.glob(os.path.join(self.included_presets_dir, self.static_files_dir_name, '*'))
            static_files_dir_in_db = os.path.join(instance_config['assnake_db'], self.dir_in_database, self.static_files_dir_name)
            os.makedirs(static_files_dir_in_db, exist_ok=True)
            # Prepare for parameter copying
            for static_file_loc in static_files_locs:
                dest_loc = os.path.join(static_files_dir_in_db, os.path.basename(static_file_loc))
                shutil.copyfile(static_file_loc, dest_loc)
            print('SUCCESSFULLY DEPLOYED PRESET PARAMETERS TO DATABASE')
        else:
            print('NOT PROPERLY CONFIGURED\nRun assnake config init')

    def cli_from_schema(self):
        pass

    def gen_click_option(self):

        instance_config = read_assnake_instance_config()
        if instance_config is not None:
            presets_glob = os.path.join(instance_config['assnake_db'], self.dir_in_database, '*.' + self.preset_file_format)
            presets = [p.split('/')[-1].replace('.json', '')
                       for p in glob.glob(presets_glob)]

            if len(presets) > 0:
                help_m = 'Preset to use. Available presets: ' + str([p.split('.')[0] for p in presets])
                default = presets[0]
            else:
                help_m = 'No presets in database!'
                default = 'No presets in database!'

            return [click.option('--preset',
                        help=help_m,
                        required=False,
                        default = default)]
        
        else: 
            return []
