import argparse
import os.path
import os
import sys
import glob
import fnmatch
from shutil import copy2, rmtree
from assnake import utils
import traceback
import parse
from assnake.api.loaders import load_df_from_db
import pandas as pd



def get_samples_from_dir(directory_with_reads, modify_name = None):
    samples_list = []  # to write samples in
    
    ext = '.fastq.gz'  # extention
    ending_variants = [{'name': 'normal', 'strands': {'R1': '_R1', 'R2': '_R2'}},
                    # {'name': 'ILLUMINA_WITH_LANE', 'strands': {'R1': '_L001_R1_001', 'R2': '_L001_R2_001'}},
                    {'name': 'ILLUMINA_001', 'strands': {'R1': '_R1_001', 'R2': '_R2_001'}},
                    {'name': 'SRA', 'strands': {'R1': '_1', 'R2': '_2'}}]  # possible endigngs of files


    for ending_variant in ending_variants:

        variant_w_ext = ending_variant['strands']['R1'] + ext
        globbing_path = os.path.join(directory_with_reads, '*' + variant_w_ext)
        files_mathching_pattern = glob.glob(globbing_path)

        # get sample names
        samples = [
            filename.split('/')[-1].replace(variant_w_ext, '') # Take last part of path (basename) and replace _R1/fastq.gz with nothing to get sample name
            for filename in files_mathching_pattern
        ]
        
        # prepare sample dicts
        samples_list += [
            {
                'name_in_run': s,
                'modified_name': modify_name(s) if modify_name is not None else s,
                
                'ending_variant_id': ending_variant['name'],
                'ending_variant_R1': ending_variant['strands']['R1'],
                'ending_variant_R2': ending_variant['strands']['R2'],
                'directory': directory_with_reads,
                'extension': ext
            } 
            for s in samples
        ]
        

    return pd.DataFrame(samples_list)


def create_links(import_dir, samples, hard = False, create_dir_if_not_exist = False, rename = False, just_print = False):
    """
    This method creates links or hard copies reads files from one directory to another. 

    :param dir_with_reads: folder where to put reads
    :param original_dir: folders that contains original reads files 
    :param sample: dictionary from get_sample_dict_from_dir method
    :param hard: if hard copying is needed or symbolic is sufficient  (False)

    :returns:

    """
    orig_wc     = '{dir_w_reads}/{name_in_run}{ending_variant}{extension}'
    new_file_wc = '{import_dir}/{name_in_dataset}{strand}{extension}'

    # make dirs
    if not os.path.isdir(import_dir):
        if create_dir_if_not_exist:
            os.makedirs(import_dir)
        else:
            print("Such dir doesn't exist and create_dir_if_not_exist is set to False")

    for sample in samples.to_dict(orient = 'records'):
        
        src_r1 = orig_wc.format(
            dir_w_reads = sample['directory'], 
            name_in_run = sample['name_in_run'],
            ending_variant = sample['ending_variant_R1'],
            extension = sample['extension']
        )
        dst_r1 = new_file_wc.format(
            import_dir = import_dir, 
            name_in_dataset = sample['fs_name'],
            strand = '_R1',
            extension = sample['extension']
        )

        src_r2 = orig_wc.format(
            dir_w_reads = sample['directory'], 
            name_in_run = sample['name_in_run'],
            ending_variant = sample['ending_variant_R2'],
            extension = sample['extension']
        )
        dst_r2 = new_file_wc.format(
            import_dir = import_dir, 
            name_in_dataset = sample['fs_name'],
            strand = '_R2',
            extension = sample['extension']
        )

        if rename:
            os.rename(src_r1, dst_r1)
            os.rename(src_r2, dst_r2)
            return

        if hard:
            copy2(src_r1, dst_r1)
            copy2(src_r2, dst_r2)
            return
        # print('---')
        # print(src_r1, dst_r1)
        # print(src_r2, dst_r2)
        os.symlink(src_r1, dst_r1)
        os.symlink(src_r2, dst_r2)

def delete_ds(dataset):
    """
    Remove assnake dataset from database
    """
    try:
        os.remove(
            '{config}/datasets/{df}/df_info.yaml'.format(config=utils.load_config_file()['assnake_db'], df=dataset))
        return (True,)
    except Exception as e:
        return (False, traceback.format_exc())

