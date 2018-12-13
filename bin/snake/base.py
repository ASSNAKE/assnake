import glob
import os
import sqlite3
import datetime


db_loc = '/data6/bio/TFM/asshole/db.sqlite3'

def save_to_db(task_id, rule_name, in_list, out_list, status ):
    
    save_str_wc = "INSERT INTO explorer_snakeruleresult VALUES (null, '{date_time}', '{task_id}', '{rule_name}', '{in_list}','{out_list}', '{status}');"
    save_str = save_str_wc.format(date_time=datetime.datetime.now(),
task_id=task_id, rule_name=rule_name, in_list=in_list, out_list=out_list, status=status)
    
    print(save_str)
    
    conn = sqlite3.connect(db_loc)
    c = conn.cursor()
    try:
        c.execute(save_str)
        conn.commit()
    except:
        print('error')
    conn.close()

configfile: 'config.yml'
    
snakefiles = os.path.join(config["software"]["snakemake_folder"], "bin/snake/")
include: snakefiles + "bowtie2.py"
include: snakefiles + "anvio.py"
include: snakefiles + "prokka.py"
include: snakefiles + "bwa.py"
include: snakefiles + "megahit.py"
include: snakefiles + "general.py"
include: snakefiles + "preprocess.py"
include: snakefiles + "marker"
include: snakefiles + "taxa.py"
include: snakefiles + "strain_finder.py"
include: snakefiles + "download.py"
include: snakefiles + "test.py"
    
