import argparse
import os.path
import os
import sys
import glob
import fnmatch
from assnake.api.fs_helpers import create_links

"Script for importing files to the assnake. Creates {prefix}/{df}/reads/{preproc}/{sample_name}"


def find_files(base, pattern):
    """
    Return list of files matching pattern in base folder.
    """
    return [n for n in fnmatch.filter(os.listdir(base), pattern) if
            os.path.isfile(os.path.join(base, n))]


def get_sample_dict_from_dir(loc, sample, variant, ext):
    temp_samples_dict = {'sample_name': sample,
                         'files': {'R1': '', 'R2': '', 'S': []},
                         'renamed_files': {'R1': '', 'R2': '', 'S': []}}

    for strand in ['R1', 'R2']:
        st = variant['strands'][strand]
        st_file = find_files(loc, sample + st + ext)
        if len(st_file) == 1:
            stripped = st_file[0].replace(variant['strands'][strand] + ext, '')
            stripped += '_' + strand + ext
            temp_samples_dict['renamed_files'][strand] = stripped
            temp_samples_dict['files'][strand] = st_file[0]

    return temp_samples_dict


def get_samples_from_dir(loc):
    """
    Searches for samples in loc. Sample should contain R1 and R2
    :param loc: location on filesystem where we should search
    :return: Returns sample dict in loc
    """
    samples_list = []
    ext = '.fastq.gz'
    end_variants = [{'name': 'normal', 'strands': {'R1': '_R1', 'R2': '_R2'}},
                    {'name': 'ILLUMINA_long', 'strands': {'R1': '_L001_R1_001', 'R2': '_L001_R2_001'}},
                    {'name': 'ILLUMINA', 'strands': {'R1': '_R1_001', 'R2': '_R2_001'}},
                    {'name': 'SRA', 'strands': {'R1': '_1', 'R2': '_2'}}]

    # Fool check
    if loc[-1] != '/':
        loc += '/'

    for variant in end_variants:
        R1 = variant['strands']['R1']  # Get first strand. For paired end data it doesn't matter.
        samples = [
            item.split("/")[-1].split(ext)[0].replace(R1, '')
            for item in glob.glob(loc + '*' + R1 + ext)
        ]
        if len(samples) > 0:
            for sample in samples:
                samples_list.append(get_sample_dict_from_dir(loc, sample, variant, ext))
            break

    return samples_list


# synonymous function was found in api.fs_helpers
'''def create_links(dir_with_reads, original_dir, sample):
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
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original_dir', type=str, help='Directory with raw reads')
    parser.add_argument('--prefix', type=str, help='Prefix on file system for dataset')
    parser.add_argument('--df', type=str, help='Name of yor dataset')
    parser.add_argument('--preproc', type=str, default='raw', help='First preprocessing')
    parser.add_argument('--dry', default=False, action='store_true',
                        help='Only print what is going to happen, without actually doing anything')
    parser.add_argument('--samples', nargs='+', help='List of sample names to import', default=[], required=False)

    args = parser.parse_args()

    # FIND SAMPLES IN DIRECTORY
    samples = get_samples_from_dir(args.original_dir)
    print(samples)
    # CREATE DIRECTORY FOR DF
    dir_with_reads = os.path.join(args.prefix, args.df, 'reads', args.preproc)
    print(dir_with_reads)
    # check if it doesn't exist
    if args.dry:
        print('Creating', dir_with_reads)
    else:
        if not os.path.isdir(dir_with_reads):
            os.makedirs(dir_with_reads)

    orig_wc = '{orig_dir}/{sample_file}'
    new_dir_wc = '{reads_dir}/{sample_name}'
    new_file_wc = '{new_dir}/{sample_file}'

    if len(args.samples) == 0:
        args.samples = [s['sample_name'] for s in samples]
    # CREATE SYMLINKS
    for s in samples:
        if len(args.samples) > 0:
            if s['sample_name'] in args.samples:

                new_dir = new_dir_wc.format(reads_dir=dir_with_reads, sample_name=s['sample_name'])

                src_r1 = orig_wc.format(orig_dir=args.original_dir, sample_file=s['files']['R1'])
                dst_r1 = new_file_wc.format(new_dir=new_dir, sample_file=s['renamed_files']['R1'])

                src_r2 = orig_wc.format(orig_dir=args.original_dir, sample_file=s['files']['R2'])
                dst_r2 = new_file_wc.format(new_dir=new_dir, sample_file=s['renamed_files']['R2'])

                if args.dry:
                    print('########', s['sample_name'], '########')
                    print('  ', new_dir)
                    print('  ', src_r1, '--->', dst_r1)
                    print('  ', src_r2, '--->', dst_r2)
                    print(s)

                else:
                    if not os.path.isdir(new_dir):
                        os.makedirs(new_dir)
                    os.symlink(src_r1, dst_r1)
                    os.symlink(src_r2, dst_r2)


if __name__ == "__main__":
    main()
