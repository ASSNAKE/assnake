from tempfile import NamedTemporaryFile
import shutil
import csv
import os
from assnake import utils
from assnake.api import fs_helpers, loaders
fields = ['sample', 'sequencing_run']


def update_fs_samples_csv(dataset, samples_to_add, preproc):
    change_all = False
    global fields
    filename = '{config}/datasets/{df}/fs_samples.tsv'.format(config=utils.load_config_file()['assnake_db'], df=dataset)
    if not os.path.exists(filename):
        df_info = loaders.load_df_from_db(dataset)
        dicts = fs_helpers.get_samples_from_dir('{}/{}/reads/raw'.format(df_info['fs_prefix'], df_info['df']))
        samples = {d['sample_name'] for d in dicts}
        with open(filename, 'w+') as file_create:
            writer = csv.DictWriter(file_create, fieldnames=fields, delimiter='\t')
            for s in samples:
                writer.writerow({'sample': s, 'sequencing_run': ''})
    if len(samples_to_add) == 0:
        change_all = True

    tempfile = NamedTemporaryFile(mode='w', delete=False)

    with open(filename, 'r') as csvfile, tempfile:
        reader = csv.DictReader(csvfile, fieldnames=fields, delimiter='\t')
        writer = csv.DictWriter(tempfile, fieldnames=fields, delimiter='\t')
        for row in reader:
            if (row['sample'] in samples_to_add) or change_all:
                row['sequencing_run'] += ' ' + preproc
            row = {'sample': row['sample'], 'sequencing_run': row['sequencing_run']}
            writer.writerow(row)
    shutil.move(tempfile.name, filename)
