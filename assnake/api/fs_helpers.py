import argparse
import os.path
import os
import sys
import glob
import fnmatch
import parse

def find_files(base, pattern):
    """
    Return list of files matching pattern in base folder.
    """
    return [n for n in fnmatch.filter(os.listdir(base), pattern) if
            os.path.isfile(os.path.join(base, n))]


def get_sample_dict_from_dir(loc, sample, variant, ext, modify_name = None):
    if modify_name is not None:
        sample_name = modify_name(sample)
    else:
        sample_name = sample

    temp_samples_dict = {'sample_name': sample_name,
                         'files': {'R1': '', 'R2': '', 'S': []},
                         'renamed_files': {'R1': '', 'R2': '', 'S': []}}


    for strand in ['R1', 'R2']:
        st = variant['strands'][strand]
        st_file = find_files(loc, sample + st + ext)
        if len(st_file) == 1:
            stripped = st_file[0].replace(variant['strands'][strand] + ext, '')
            stripped = sample_name + '_' + strand + ext
            temp_samples_dict['renamed_files'][strand] = stripped
            temp_samples_dict['files'][strand] = st_file[0]

    return temp_samples_dict


def get_samples_from_dir(loc, modify_name = None):
    """
    Searches for samples in loc. Sample should contain R1 and R2
    :param loc: location on filesystem where we should search
    :return: Returns list of sample dicts in loc
    """
    samples_list = []
    ext = '.fastq.gz'
    end_variants = [{'name': 'normal', 'strands': {'R1': '_R1', 'R2': '_R2'}},
                    {'name': 'ILLUMINA_1', 'strands': {'R1': '_L001_R1_001', 'R2': '_L001_R2_001'}},
                    {'name': 'ILLUMINA', 'strands': {'R1': '_R1_001', 'R2': '_R2_001'}},
                    {'name': 'SRA', 'strands': {'R1': '_1', 'R2': '_2'}}]

    # Fool check
    if loc[-1] != '/':
        loc += '/'

    for variant in end_variants:
        R1 = variant['strands']['R1']  # Get just the first strand. For paired end data it doesn't matter.
        samples = [
            item.split("/")[-1].split(ext)[0].replace(R1, '')
            for item in glob.glob(loc + '*' + R1 + ext)
        ]
        if len(samples) > 0:
            for sample in samples:
                samples_list.append(get_sample_dict_from_dir(loc, sample, variant, ext, modify_name))
            break

    return samples_list

def create_links(dir_with_reads, original_dir, sample):
    orig_wc = '{orig_dir}/{sample_file}'
    new_dir_wc = '{reads_dir}/{sample_name}'
    new_file_wc = '{new_dir}/{sample_file}'

    new_dir = new_dir_wc.format(reads_dir = dir_with_reads, sample_name = sample['sample_name'])
    if not os.path.isdir(new_dir):
        os.makedirs(new_dir)
        
    src_r1 = orig_wc.format(orig_dir = original_dir, sample_file = sample['files']['R1'])
    dst_r1 = new_file_wc.format(new_dir=new_dir, sample_file=sample['renamed_files']['R1'])

    src_r2 = orig_wc.format(orig_dir = original_dir, sample_file = sample['files']['R2'])
    dst_r2 = new_file_wc.format(new_dir=new_dir, sample_file=sample['renamed_files']['R2'])

    os.symlink(src_r1, dst_r1)
    os.symlink(src_r2, dst_r2)
    # try:
    #     os.symlink(src_r1, dst_r1)
    #     os.symlink(src_r2, dst_r2)
    # except:
    #     print(sample, 'ERROR')