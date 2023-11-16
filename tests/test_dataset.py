import pytest
from assnake.utils.fs_helpers import get_samples_from_dir, create_links
from tests.util_for_test import random_path, random_file_name
import random
import pathlib
import os


dest2check = [random_path() for _ in range(1)]
smpls_dict = dict()
@pytest.mark.smoke
def test_dataset_finding(test_assn_data):
    '''
        testing get_samples_from_dir
    '''
    global smpls_dict
    buf_dir= test_assn_data['written_files']
    source_reads2check, file_abs_paths, samples_names = buf_dir['source_reads2check'], buf_dir['file_abs_paths'], buf_dir['samples_names']
    path_to_smpls_dir = test_assn_data['smpls']
    smpls = list()
    os_find = list()
    for loc in source_reads2check:
        respond = get_samples_from_dir(str(path_to_smpls_dir / loc))
        smpls_dict.update({loc: respond})
        #print('respond: ',respond)
        smpls +=[i['sample_name'] for i in respond]
        #print(get_samples_from_dir(str(path_to_smpls_dir / loc)))
        #print([i['sample_name'] for i in respond])
        #print(smpls)
        #print(samples_names, len(set(samples_names)))
        os_find += os.listdir(str(path_to_smpls_dir / loc))
        # print(os_find)
        # print(smpls)
    assert set(smpls) == set(samples_names)

@pytest.mark.smoke
@pytest.mark.dataset_api
@pytest.mark.parametrize('dir_with_reads', dest2check)
def test_dataset_import(test_assn_data, dir_with_reads):
    global smpls_dict
    buf_dir = test_assn_data['written_files']
    source_reads2check, file_abs_paths, samples_names = buf_dir['source_reads2check'], buf_dir['file_abs_paths'], \
                                                       buf_dir['samples_names']
    path_to_smpls_dir = test_assn_data['smpls']
    os.makedirs(str(path_to_smpls_dir / dir_with_reads))
    renames2check = list()
    for i, original_reads in enumerate(source_reads2check):
        for dict in get_samples_from_dir(str(path_to_smpls_dir / original_reads)):
            create_links(str(path_to_smpls_dir / dir_with_reads), original_reads, dict)
            renames2check+=['{sample_file}'.format(sample_file=dict['renamed_files']['R2'])]
            renames2check += ['{sample_file}'.format(sample_file=dict['renamed_files']['R1'])]
    assert set(os.listdir(str(path_to_smpls_dir / dir_with_reads))) == set(renames2check)



