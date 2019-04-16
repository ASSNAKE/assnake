# This is a small test suite you should run in order to test your installation of ASSNAKE
import os
# import glob
import json
import pandas as pd

# Configuration
assnake_install_dir = '/data6/bio/TFM/pipeline/assnake/'
configfile: os.path.join(assnake_install_dir, 'config.yml')
include: os.path.join(assnake_install_dir, 'bin/snake/base.py')

# Load info about results
results = glob(os.path.join(assnake_install_dir, './results/*/*.json'))
results_dicts = []
for res in results:
    with open(res) as result_file:
        results_dicts.append(json.load(result_file))

results = pd.DataFrame(results_dicts)

print(results)

test_df_prefix = config['test_df_prefix']
test_df_name = config['test_df_name']

samples_wc = '{prefix}/{df}/reads/{preproc}/{sample}'
samples = [r.split('/')[-1] for r in glob(samples_wc.format(
    prefix = test_df_prefix,
    df = test_df_name,
    preproc = 'raw',
    sample = '*'
))]

count_wc = results.loc[results['short_name'] == 'count']['out_str_wc'].values[0]
count_wc = count_wc.replace('strand', '{strand}')
print(count_wc)

rule test_count: 
    input: expand(count_wc.format(prefix = test_df_prefix, df = test_df_name, preproc='raw', sample = samples[0]), strand = ['R1', 'R2'])

tmtic_wc = results.loc[results['short_name'] == 'tmtic']['out_str_wc'].values[0]
rule test_tmtic: 
    input: tmtic_wc.format(prefix = test_df_prefix, params = 'def1', df = test_df_name, preproc='raw', sample = samples[1])
