import argparse
import os.path
import os
import sys
import glob
import fnmatch
import parse
from assnake.api.loaders import load_dfs_from_db

def find_files(base, pattern):
    """
    Return list of files matching pattern in base folder.
    """
    return [n for n in fnmatch.filter(os.listdir(base), pattern) if
            os.path.isfile(os.path.join(base, n))]


def get_sample_dict_from_dir(loc, sample_name, variant, ext, modify_name=lambda arg: arg):
    '''

    :param loc:
    :param sample_name:
    :param variant:
    :param ext:
    :param modify_name:
    :return:
    '''
    # DONE в 1 строчку

    sample_name = modify_name(sample_name)
    temp_samples_dict = {'sample_name': sample_name,
                         'files': {'R1': '', 'R2': '', 'S': []},
                         'renamed_files': {'R1': '', 'R2': '', 'S': []}
                         }

    for strand in ['R1', 'R2']:
        st = variant['strands'][strand]
        st_file = find_files(loc, sample_name + st + ext)
        if len(st_file) == 1:
            # DONE переписать через format
            stripped = st_file[0].replace(variant['strands'][strand] + ext, '')
            stripped += '{sample_name}_{stard}{ext}'.format(sample_name=sample_name, strand=strand, ext=ext)
            temp_samples_dict['renamed_files'][strand] = stripped
            temp_samples_dict['files'][strand] = st_file[0]

    return temp_samples_dict


def get_samples_from_dir(loc, modify_name=None):
    """
    Searches for samples in loc. Sample should contain R1 and R2
    :param loc: location on filesystem where we should search
    :return: Returns list of sample dicts in loc
    """
    samples_list = []  # to write samples in
    ext = '.fastq.gz'  # extention
    end_variants = [{'name': 'normal', 'strands': {'R1': '_R1', 'R2': '_R2'}},
                    # {'name': 'ILLUMINA_1', 'strands': {'R1': '_L001_R1_001', 'R2': '_L001_R2_001'}},
                    {'name': 'ILLUMINA', 'strands': {'R1': '_R1_001', 'R2': '_R2_001'}},
                    {'name': 'SRA', 'strands': {'R1': '_1', 'R2': '_2'}}]  # possible endigngs of files

    # Fool check
    if loc[-1] != '/':
        loc += '/'

    for variant in end_variants:
        R1 = variant['strands']['R1']  # Get just the first strand. For paired end data it doesn't matter.
        samples = [
            item[item.rfind('/')+1:item.rfind(R1+ext)]
            for item in glob.glob(loc + '*' + R1 + ext)
        ]
        if len(samples) > 0:
            for sample in samples:
                samples_list.append(get_sample_dict_from_dir(loc, sample, variant, ext, modify_name))
            break

    return samples_list


# TODO везде документацию,
# DONE все класть {fs_prefix}/{df}/reads/{preproc}/{sample}_{strand}.fastq.gz
def create_links(dir_with_reads, original_dir, sample):
    """
    :param dir_with_reads:  куда класть
    :param original_dir: откуда
    :param sample: dictionary от get samples dict from dir
    :return:

    """

    orig_wc = '{orig_dir}/{sample_file}'
    new_file_wc = '{new_dir}/{sample_name}_{sample_file}'

    if not os.path.isdir(dir_with_reads):
        os.makedirs(dir_with_reads)

    src_r1 = orig_wc.format(orig_dir=original_dir, sample_file=sample['files']['R1'])
    dst_r1 = new_file_wc.format(new_dir=dir_with_reads, sample_name = sample['sample_name'],  sample_file=sample['renamed_files']['R1'])

    src_r2 = orig_wc.format(orig_dir=original_dir, sample_file=sample['files']['R2'])
    dst_r2 = new_file_wc.format(new_dir=dir_with_reads, sample_name = sample['sample_name'],  sample_file=sample['renamed_files']['R2'])

    os.symlink(src_r1, dst_r1)
    os.symlink(src_r2, dst_r2)
    # try:
    #     os.symlink(src_r1, dst_r1)
    #     os.symlink(src_r2, dst_r2)
    # except:
    #     print(sample, 'ERROR')
