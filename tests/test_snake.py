# This is a small test suite you should run in order to test your installation of ASSNAKE
import os
import json
import pandas as pd

# Configuration
include: os.path.join(config['assnake_install_dir'], 'bin/snake/base.py')


# Load info about results
results = glob(os.path.join(config['assnake_install_dir'], './results/*/*.json'))
results_dicts = []
for res in results:
    with open(res, encoding="utf-8") as result_file:
        results_dicts.append(json.load(result_file))
results = pd.DataFrame(results_dicts)


test_df_prefix = config['test_df_prefix']
test_df_name = config['test_df_name']

samples_wc = '{prefix}/{df}/reads/{preproc}/{sample}'
samples = [r.split('/')[-1] for r in glob(samples_wc.format(
    prefix = test_df_prefix,
    df = test_df_name,
    preproc = 'raw',
    sample = '*'
))]

samp_ind = 3
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


i_want = count_locs + [tmtic_req]

rule main_test:
    input: i_want

