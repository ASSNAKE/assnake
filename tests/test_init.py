import pytest
from assnake.api.init_config import fill_and_write_config
import pathlib
from util_for_test import random_path


args2check = [[random_path() for _ in range(6)] for __ in range(10)]



@pytest.mark.init
@pytest.mark.smoke
@pytest.mark.parametrize('paths2check', args2check)
def test_api_init_start_simple_random(test_assn_data, paths2check):
    path_to_db = test_assn_data['db']
    paths2compare = list(map(lambda rel: '{pre}/{rel}'.format(pre=path_to_db, rel=rel), paths2check))
    print(paths2compare)
    fill_and_write_config(*paths2compare)
    paths2compare[5] = paths2compare[5]+'/config.yaml'
    for i in range(1,6):
        assert pathlib.Path(paths2compare[i]).exists()


@pytest.mark.init
def test_api_init_start_nested(test_assn_data):
    path_to_db = test_assn_data['db']
    paths2check = [
        './fictive/../target',
        './fictive/../target/1',
        './fictive/../target/2',
        './fictive/../target/3',
        './fictive/../target/4',
        './fictive/../target/5'
    ]
    paths2compare = list(map(lambda rel: '{pre}/{rel}'.format(pre=path_to_db, rel=rel), paths2check))
    fill_and_write_config(*paths2compare)

    assert (path_to_db / 'target').exists()
    for i in range(1,5):
        assert (path_to_db / ('target/'+str(i))).exists()
    assert (path_to_db / 'target/5/config.yaml').exists()