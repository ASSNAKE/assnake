# This is a small test suite you should run in order to test your installation of ASSNAKE
import os
import json
import pandas as pd
import yaml

# LOAD WC CONFIG
wc_config = None
wc_config_loc = os.path.join(dir_of_this_file, '../wc_config.yaml')
with open(wc_config_loc, 'r') as stream:
    try:
        wc_config = yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Configuration
include: os.path.join(config['assnake_install_dir'], 'bin/snake/base.py')


# Load info about results
results = glob(os.path.join(config['assnake_install_dir'], './results/*/*.json'))
results_dicts = []
for res in results:
    with open(res, encoding="utf-8") as result_file:
        results_dicts.append(json.load(result_file))
results = pd.DataFrame(results_dicts)

print('installation_dir', config['assnake_install_dir'])
test_df_prefix = config['test_df_prefix']
test_df_name = config['test_df_name']

samples_wc = '{prefix}/{df}/reads/{preproc}/{sample}'
samples = [r.split('/')[-1] for r in glob(samples_wc.format(
    prefix = test_df_prefix,
    df = test_df_name,
    preproc = 'raw',
    sample = '*'
))]

samp_ind = 4
count_wc = results.loc[results['short_name'] == 'count']['out_str_wc'].values[0]
count_wc = count_wc.replace('strand', '{strand}')

count_locs = expand(count_wc.format(prefix = test_df_prefix, 
                df = test_df_name, preproc='raw', 
                sample = samples[samp_ind]), strand = ['R1', 'R2'])

tmtic_wc = results.loc[results['short_name'] == 'tmtic']['out_str_wc'].values[0]
tmtic_req = tmtic_wc.format(prefix = test_df_prefix, 
                            params = 'def1', 
                            df = test_df_name, 
                            preproc='raw', 
                            sample = samples[samp_ind])

mp2_wc = results.loc[results['short_name'] == 'metaphlan2']['out_str_wc'].values[0]
centr_wc = results.loc[results['short_name'] == 'centrifuge']['out_str_wc'].values[0]

mp2_loc = mp2_wc.format(prefix = test_df_prefix, 
                            df = test_df_name, 
                            preproc='raw', 
                            sample = samples[samp_ind])
centr_loc = centr_wc.format(prefix = test_df_prefix, 
                            params = 'def',
                            df = test_df_name, 
                            preproc='raw', 
                            sample = samples[samp_ind])

i_want = count_locs + [tmtic_req, mp2_loc, centr_loc]




rule main_test:
    input: i_want

