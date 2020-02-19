"""
here we make configurations for testing like mocks, stubs, fixtures -- I don't know, frankly speaking, but it have to help in testing our API
"""
import random, os
from pathlib import Path
# from typing import Dict

from tests.util_for_test import random_file_name, random_path
import pytest
import shutil
file_abs_paths = list()
source_reads2check = list()
samples_names = list()
@pytest.fixture(scope="session", autouse=True)
def test_assn_data(tmp_path_factory):
    '''
    It`s built-in fixture of pytest module, which produces temporary directory, which could be accessed in test scripts
    as def test(test_assn_data): .. return;
    '''
    assn_db_testf = tmp_path_factory.mktemp('test_assn_data')
    assn_samples_dir = tmp_path_factory.mktemp('test_samples_dir')
    dirs = {'db': assn_db_testf, 'smpls': assn_samples_dir}

    endings = ['_R{strand}.fastq.gz', '_R{strand}_001.fastq.gz', '_{strand}.fastq.gz']

    global source_reads2check
    source_reads2check = [random_path() for _ in range(10)]

    global samples_names
    file_rel_paths = list()
    for n in range(10):
        os.makedirs(assn_samples_dir / source_reads2check[n])
        for k in range(n):
            buff = random_file_name()
            samples_names += [buff]
            pre_file_name = buff + random.choice(endings)
            pre_file_name = [pre_file_name.format(strand=1), pre_file_name.format(strand=2)]
            file_rel_paths += [os.path.join(source_reads2check[n], pre_file_name[0]),
                               os.path.join(source_reads2check[n], pre_file_name[1])]
    global file_abs_paths
    file_abs_paths = list(map(lambda rel: '{pre}/{rel}'.format(pre=assn_samples_dir, rel=rel), file_rel_paths))

    for i in file_abs_paths:
        with open(i, "w+") as f:
            f.write('1')
    dirs.update({'written_files':{'source_reads2check':source_reads2check, 'file_abs_paths':file_abs_paths, 'samples_names': samples_names}})
    yield dirs
    shutil.rmtree(dirs['db'])
    shutil.rmtree(dirs['smpls'])



