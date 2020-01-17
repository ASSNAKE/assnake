import glob
import os
import sqlite3
import datetime
import yaml

wc_config_loc = os.path.join(config['assnake_install_dir'], 'wc_config.yaml')
wc_config = {}
with open(wc_config_loc, 'r') as stream:
    try:
        wc_config = yaml.load(stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        print(exc)

db_loc = '/data6/bio/TFM/asshole/db.sqlite3'

def save_to_db(task_id, rule_name, in_list, out_list, status ):
    
    save_str_wc = "INSERT INTO explorer_snakeruleresult VALUES (null, '{date_time}', '{task_id}', '{rule_name}', '{in_list}','{out_list}', '{status}');"
    save_str = save_str_wc.format(date_time=datetime.datetime.now(),
task_id=task_id, rule_name=rule_name, in_list=in_list, out_list=out_list, status=status)
    
    print(save_str)
    
    conn = sqlite3.connect(db_loc)
    c = conn.cursor()
    c.execute(save_str)
    conn.commit()
    conn.close()
    
    
snakefiles = './'
modules    = '../../../snake/modules/'

# == MAPPING ==
include: modules + 'bwa/bwa.py'
include: modules + 'bbmap/bbmap.py'
include: modules + 'coverage/profile.py'

# == QUALITY CONTROL ==
include: modules + 'count/count.py'
include: modules + 'fastqc/fastqc.py'
include: modules + 'preprocessing/gz_operations.py'
include: modules + 'remove_human_bbmap/remove_human_bbmap.py'
include: modules + 'trim_galore/trim_galore.py'
include: modules + 'minlen/minlen_bbmap.py'
include: modules + 'trimmomatic/trimmomatic.py'
include: modules + 'multiqc/multiqc.py'

# == TAXONOMY ==
include: modules + 'metaphlan2/metaphlan2.py'
include: modules + 'centrifuge/centrifuge.py'

include: modules + 'anvio/anvio.py'
include: modules + 'blastn/blastn.py'
include: modules + 'metawrap_classify_bins/metawrap_classify_bins.py'
include: modules + 'checkm/checkm.py'
include: modules + 'strain_finder/strain_finder.py'
include: modules + 'cat_bat/cat_bat.py'
include: modules + 'kraken2/kraken2.py'

# == ASSEMBLY ==
include: modules + 'metaspades/metaspades.py'
include: modules + 'megahit/megahit_cross.py'

# == BINNING ==
include: modules + 'maxbin2/maxbin2.py'
include: modules + 'metabat2/metabat2.py'

# == 16S ==
include: modules + 'dada2/dada2.py'
include: modules + 'picrust2/picrust2.py'
# include: modules + 'mothur/mothur.py'







include: snakefiles + "bins.py"


include: snakefiles + "humann2.py"
include: snakefiles + "qiime2.py"
include: snakefiles + "fasta_operations.py"


include: snakefiles + "general.py"
include: snakefiles + "preprocess.py"
include: snakefiles + "taxa.py"
include: snakefiles + "download.py"
include: snakefiles + "find_fungi.py"
    
