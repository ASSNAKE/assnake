configfile: '/data4/bio/fedorov/assnake_dev/config.yml'
include: "/data4/bio/fedorov/assnake_dev/bin/snake/base.py"

import glob

samples = [s.split('/')[-1] for s in glob.glob('/data6/bio/TFM/datasets/mg_test/reads/raw/*')]


i_want = expand('{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped.sam', 
        prefix = '/data6/bio/TFM/datasets',
        df = 'mg_test',
        params='def1',
        path = 'db/MEGARes/v1.01',
        seq_set_id = 'megares',
        preproc = 'raw',
        sample = samples)
rule map_on_smth:
    input: i_want
