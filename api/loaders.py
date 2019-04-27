import os
import glob

import pandas as pd
import numpy as np
import bb_stats

import yaml

data_dir = "/data6/bio/TFM/pipeline/"
DF_DIR = data_dir + 'datasets/{df}/'
SAMPLE_DIR = DF_DIR + 'reads/{preproc}/{sample}'
FQGZ_LOC = SAMPLE_DIR +  '/{sample}_{strand}.fastq.gz'
COUNT_LOC = SAMPLE_DIR + '/profile/{sample}_{strand}.count'

SAMPLES_META_LOC = DF_DIR + '{df}_samples.tsv'
SOURCES_META_LOC = DF_DIR + '{df}_sources.tsv'

pipeline = '/data6/bio/TFM/pipeline/'


DF_DIR_PREFIX = '{prefix}/{df}/'
SAMPLE_DIR_PREFIX = DF_DIR_PREFIX + 'reads/{preproc}/{sample}'
FQGZ_LOC_PREFIX = SAMPLE_DIR_PREFIX +  '/{sample}_{strand}.fastq.gz'
COUNT_LOC_PREFIX = SAMPLE_DIR_PREFIX + '/profile/{sample}_{strand}.count'

SAMPLES_META_LOC = DF_DIR_PREFIX + '{df}_samples.tsv'
SOURCES_META_LOC = DF_DIR_PREFIX + '{df}_sources.tsv'

## {{{ http://code.activestate.com/recipes/578019/ (r15)

"""
Bytes-to-human / human-to-bytes converter.
Based on: http://goo.gl/kTQMs
Working with Python 2.x and 3.x.

Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
License: MIT
"""

# see: http://goo.gl/kTQMs
SYMBOLS = {
    'customary'     : ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
    'customary_ext' : ('byte', 'kilo', 'mega', 'giga', 'tera', 'peta', 'exa',
                       'zetta', 'iotta'),
    'iec'           : ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
    'iec_ext'       : ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi',
                       'zebi', 'yobi'),
}

def bytes2human(n, format='%(value).1f %(symbol)s', symbols='customary'):
    """
    Convert n bytes into a human readable string based on format.
    symbols can be either "customary", "customary_ext", "iec" or "iec_ext",
    see: http://goo.gl/kTQMs

      >>> bytes2human(0)
      '0.0 B'
      >>> bytes2human(0.9)
      '0.0 B'
      >>> bytes2human(1)
      '1.0 B'
      >>> bytes2human(1.9)
      '1.0 B'
      >>> bytes2human(1024)
      '1.0 K'
      >>> bytes2human(1048576)
      '1.0 M'
      >>> bytes2human(1099511627776127398123789121)
      '909.5 Y'

      >>> bytes2human(9856, symbols="customary")
      '9.6 K'
      >>> bytes2human(9856, symbols="customary_ext")
      '9.6 kilo'
      >>> bytes2human(9856, symbols="iec")
      '9.6 Ki'
      >>> bytes2human(9856, symbols="iec_ext")
      '9.6 kibi'

      >>> bytes2human(10000, "%(value).1f %(symbol)s/sec")
      '9.8 K/sec'

      >>> # precision can be adjusted by playing with %f operator
      >>> bytes2human(10000, format="%(value).5f %(symbol)s")
      '9.76562 K'
    """
    n = int(n)
    if n < 0:
        raise ValueError("n < 0")
    symbols = SYMBOLS[symbols]
    prefix = {}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i+1)*10
    for symbol in reversed(symbols[1:]):
        if n >= prefix[symbol]:
            value = float(n) / prefix[symbol]
            return format % locals()
    return format % dict(symbol=symbols[0], value=n)

def human2bytes(s):
    """
    Attempts to guess the string format based on default symbols
    set and return the corresponding bytes as an integer.
    When unable to recognize the format ValueError is raised.

      >>> human2bytes('0 B')
      0
      >>> human2bytes('1 K')
      1024
      >>> human2bytes('1 M')
      1048576
      >>> human2bytes('1 Gi')
      1073741824
      >>> human2bytes('1 tera')
      1099511627776

      >>> human2bytes('0.5kilo')
      512
      >>> human2bytes('0.1  byte')
      0
      >>> human2bytes('1 k')  # k is an alias for K
      1024
      >>> human2bytes('12 foo')
      Traceback (most recent call last):
          ...
      ValueError: can't interpret '12 foo'
    """
    init = s
    num = ""
    while s and s[0:1].isdigit() or s[0:1] == '.':
        num += s[0]
        s = s[1:]
    num = float(num)
    letter = s.strip()
    for name, sset in SYMBOLS.items():
        if letter in sset:
            break
    else:
        if letter == 'k':
            # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
            sset = SYMBOLS['customary']
            letter = letter.upper()
        else:
            raise ValueError("can't interpret %r" % init)
    prefix = {sset[0]:1}
    for i, s in enumerate(sset[1:]):
        prefix[s] = 1 << (i+1)*10
    return int(num * prefix[letter])


## end of http://code.activestate.com/recipes/578019/ }}}



def load_count(prefix, df, preproc, sample):
    strands = ['R1', 'R2']
    count_loc1 = COUNT_LOC_PREFIX.format(prefix=prefix, df=df, preproc=preproc, sample=sample, strand=strands[0])
    count_loc2 = COUNT_LOC_PREFIX.format(prefix=prefix, df=df, preproc=preproc, sample=sample, strand=strands[1])
    
    reads = -1
    bps = -1
    
    try:
        with open(count_loc1, 'r') as f:
            line = f.readline().rstrip()
            reads += int(line.split(' ')[0])
            bps += int(line.split(' ')[1])
        with open(count_loc2, 'r') as f:
            line = f.readline().rstrip()
            reads += int(line.split(' ')[0])
            bps += int(line.split(' ')[1])
    except:
        print('error loading counts: ' + sample)
        return {'reads': -1, 'bps': -1}
    
    return {'reads': reads+1, 'bps': bps+1}
        

def load_dfs_from_db(db_loc):
    """
    Returns list of dictionaries with info about datasets from fs database.
    Mandatory fields: df, prefix
    """
    dfs = []
    df_info_locs = glob.glob(db_loc+'/datasets/*/df_info.yaml')
    
    for df_info in df_info_locs:
        with open(df_info, 'r') as stream:
            try:
                dfs.append(yaml.load(stream))
            except yaml.YAMLError as exc:
                print(exc)
    return dfs


def load_sample(prefix, df, preproc, sample):
    strands = ['R1', 'R2']
    
    final_preproc = ''
    size = 0
    containers = []
    preprocs = []
    if preproc == 'longest':
        preprocs = [p.split('/')[-2] for p in 
                    glob.glob(SAMPLE_DIR_PREFIX.format(prefix=prefix, df=df, preproc='*', sample=sample))]
    else:
        preprocs = [p.split('/')[-2] for p in 
                    glob.glob(SAMPLE_DIR_PREFIX.format(prefix=prefix, df=df, preproc=preproc, sample=sample))]
    for p in preprocs:
        r1 = FQGZ_LOC_PREFIX.format(prefix=prefix, df=df, preproc=p, sample=sample, strand=strands[0])
        r2 = FQGZ_LOC_PREFIX.format(prefix=prefix, df=df, preproc=p, sample=sample, strand=strands[1])

        if os.path.isfile(r1) and os.path.isfile(r2):
            containers.append(p)
            if len(p) > len(final_preproc):
                final_preproc = p
                size = os.path.getsize(r1) + os.path.getsize(r2)

    return {'df':df, 'fs_name':sample, 'sample':sample,  'preproc':final_preproc, 
                         'preprocs':containers, 
                         'size': bytes2human(size, symbols='iec'), 'bytes': size,
                         **load_count(prefix, df, final_preproc, sample)}

def samples_in_df(df, db_loc):
    df_info = None
    df_info_loc = db_loc+'/datasets/'+df+'/df_info.yaml'
    
    with open(df_info_loc, 'r') as stream:
        try:
            df_info = (yaml.load(stream))
        except yaml.YAMLError as exc:
            print(exc)
                       
    samples = df_full_info(df_info['fs_prefix'], df_info['df'])
    return samples

def load_mg_samples_in_df_fs(db_loc, df):
    df_info = None
    df_info_loc = db_loc+'/datasets/'+df+'/df_info.yaml'
    
    with open(df_info_loc, 'r') as stream:
        try:
            df_info = (yaml.load(stream))
        except yaml.YAMLError as exc:
            print(exc)
                       
    samples = mg_samples_for_df_fs(df_info['fs_prefix'], df_info['df'])
    
    samples_info = load_mg_samples_in_df(df, db_loc, 'pandas')
    
    for sample in samples:
        biospecimen = list(samples_info.loc[samples_info['fs_name'] == sample['name_on_fs']]['biospecimen'])
        if len(biospecimen) == 1:
            sample.update({'biospecimen' : biospecimen[0]})
        else:
            print('biospecimen: ' + sample)
    
    return samples


def load_mg_files(prefix, df, preproc, sample):
    strands = ['R1', 'R2']
    files = []
    
    for strand in strands:
        count = COUNT_LOC_PREFIX.format(prefix=prefix, df=df, preproc=preproc, sample=sample, strand=strand)
        fqgz = FQGZ_LOC_PREFIX.format(prefix=prefix, df=df, preproc=preproc, sample=sample, strand=strand)
        
        reads = -1
        bps = -1
        size = 0
        try:
            with open(count, 'r') as f:
                line = f.readline().rstrip()
                reads = int(line.split(' ')[0])
                bps = int(line.split(' ')[1])
        except:
            print('error: ' + sample)
                
        files.append({'reads': reads, 'bps': bps, 'size': os.path.getsize(fqgz), 'strand': strand})
    return files

def mg_samples_for_df_fs(prefix, df):
    
    strands = ['R1', 'R2']
    df_dir = DF_DIR_PREFIX.format(df = df, prefix=prefix)
    
    if os.path.isdir(df_dir): # Fool's check
        samples = set([s.split('/')[-1] for s in glob.glob(SAMPLE_DIR_PREFIX.format(prefix=prefix,df=df,preproc='*',sample='*'))])    
        sample_dicts = []
        
        print('TOTAL SAMLES: ', len(samples))
        
        for i, s in enumerate(samples):
            sample_dict = {
                'df': df,
                'name_on_fs': s,
                'containers': []
            }
            preprocs = set([s.split('/')[-2] for s in glob.glob(SAMPLE_DIR_PREFIX.format(prefix=prefix,df=df,preproc='*',sample=s))])
            containers = []
            
            for p in preprocs:
                r1 = FQGZ_LOC_PREFIX.format(prefix=prefix, df=df, preproc=p, sample=s, strand=strands[0])
                r2 = FQGZ_LOC_PREFIX.format(prefix=prefix, df=df, preproc=p, sample=s, strand=strands[1])

                if os.path.isfile(r1) and os.path.isfile(r2):
                    containers.append({'preprocessing':p,
                                 'files':load_mg_files(prefix, df, p, s)})
                    
            sample_dict['containers'] = containers
            sample_dicts.append(sample_dict)

        return sample_dicts

def load_sources_in_df(df, db_loc, return_as='dict'):
    """
    Returns dict/dataframe with techical info about sources in dataset
    """
    loc = db_loc +'/datasets/{df}/sources.tsv'
    loc = loc.format(df=df)
    sources = pd.read_csv(loc, sep = '\t')
    if return_as == 'dict':
        return sources.to_dict(orient='records')
    elif return_as == 'pandas':
        return sources
    return 0

def load_biospecimens_in_df(df, db_loc, return_as='dict'):
    """
    Returns dict/dataframe with techical info about biospecimens in dataset
    """
    loc = db_loc +'/datasets/{df}/biospecimens.tsv'
    loc = loc.format(df=df)
    biospecimens = pd.read_csv(loc, sep = '\t')
    biospecimens = biospecimens.fillna('')
    if return_as == 'dict':
        return biospecimens.to_dict(orient='records')
    elif return_as == 'pandas':
        return biospecimens
    return 0

def load_mg_samples_in_df(df, db_loc, return_as='dict'):
    """
    Returns dict/dataframe with techical info about mg_samples in dataset
    """
    loc = db_loc +'/datasets/{df}/mg_samples.tsv'
    loc = loc.format(df=df)
    mg_samples = pd.read_csv(loc, sep = '\t')
    if return_as == 'dict':
        return mg_samples.to_dict(orient='records')
    elif return_as == 'pandas':
        return mg_samples
    return 0
 

def df_full_info(prefix, df, preprocessing = 'longest'):
    """
    Returns dict with samples [{'df': df, 'preproc': preproc, 'sample': sample}]
    
    params: 
        df - name of dataset
        preprocessing - longest/newest
    
    """
    
    strands = ['R1', 'R2']
    df_dir = DF_DIR_PREFIX.format(df = df, prefix=prefix)
    
    if os.path.isdir(df_dir): # Fool's check
        samples = set([s.split('/')[-1] for s in glob.glob(SAMPLE_DIR_PREFIX.format(prefix=prefix,df=df,preproc='*',sample='*'))])    
        sample_dicts = []
        
        for s in samples:
            sample_dicts.append(load_sample(prefix, df, preprocessing, s))
        return sample_dicts

def samples_to_pd(samples):
    meta_df = pd.DataFrame(columns=['fs_name', 'df', 'preproc', 'size', 'bytes', 'reads','bps'])
    for s in samples:
        meta_df.loc[s['sample']]=[s['fs_name'], s['df'], s['preproc'], s['size'], s['bytes'], s['reads'],s['bps']]
    return meta_df.sort_values('fs_name')

def load_samples_metadata(prefix, df):
    meta_loc = SAMPLES_META_LOC.format(prefix=prefix,df=df)
    meta = pd.read_csv(meta_loc, sep='\t')
    
    sources_loc = SOURCES_META_LOC.format(prefix=prefix,df=df)
    if os.path.isfile(sources_loc):
        sources = pd.read_csv(sources_loc, sep='\t')
        meta = meta.merge(sources, left_on='source', right_on='source' )
    return meta


def load_mp2(prefix, samples, level='s__', org='Bacteria', index_by = 'fs_name'):
    dfs = []
    loc_wc = '{prefix}/{df}/taxa/{preproc}/mp2__def/{sample}/{sample}.mp2'

    for s in samples:
        
        loc = loc_wc.format(prefix=prefix, df = s['df'], preproc=s['preproc'], sample = s['fs_name'])
        if os.path.isfile(loc):
            try:
                mp2_df = pd.read_csv(loc, sep='\t')
                mp2_df = mp2_df.set_index('#SampleID')
                mp2_df.columns = [s[index_by]]
                mp2_df.index.name = ''
                mp2_df = mp2_df.T
                dfs.append(mp2_df)
            except Exception as e:
                print('ERROR: ' + loc)
                print(e)
    mp2_combined = pd.concat(dfs, axis=0)
    mp2_combined = mp2_combined.fillna(0)

    levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 't__']
    org = 'Bacteria'

    throw_away = levels[levels.index(level)+1]
    cols = list(mp2_combined.columns)
    drop_cols = []
    for col in cols:
        if not (throw_away not in col and level in col and levels[0]+org in col):
            drop_cols.append(col)

    mp2_combined = mp2_combined.drop(drop_cols, axis=1)
    
    return mp2_combined


def load_resanal_reports(samples, level = 'Mechanism', norm = True):
    if level == 'Mechanism':
        report_wc = pipeline + 'datasets/{df}/resanal/{sample}/{preproc}/mech_fp.tsv'
    elif level == 'Group':
        report_wc = pipeline + 'datasets/{df}/resanal/{sample}/{preproc}/group_fp.tsv'

    
    dfs = []
    meta = []
    
    for s in samples:
        if os.path.isfile(report_wc.format(df = s['df'], preproc=s['preproc'], sample = s['fs_name'])):
            
            df = pd.read_csv(report_wc.format(df = s['df'], preproc=s['preproc'], sample = s['fs_name']), sep='\t')
            df = df.T
            df.columns = list(df.loc[level])
            df = df.drop('Sample')
            df = df.drop(level)
            df = df.rename(index={'Hits': s['sample']})
            
            if norm:
                count = load_sample(df = s['df'], preproc=s['preproc'], sample = s['fs_name'])
                df = df.divide(count['reads'])
                
            dfs.append(df)
        
    combined = pd.concat(dfs, axis=0)
    combined = combined.fillna(0)
    return combined

def load_hm2(prefix, samples, dbs = 'chocophlan__uniref90', index_by = 'fs_name', norm = False, modifier='unstratified'):
    dfs = []
    
    db_nucl = 'chocophlan'
    db_protein = 'uniref90'
    loc_wc = '{prefix}/{df}/humann2/{dbs}/{sample}/{preproc}/{sample}_pathabundance{modifier}.tsv'
    
    if modifier != '':
        modifier = '_' + modifier
        
    for s in samples:
        loc = loc_wc.format(prefix=prefix,
                            df = s['df'], 
                            preproc=s['preproc'], 
                            sample = s['fs_name'],
                            dbs = dbs,
                           modifier=modifier)
        if os.path.isfile(loc):
            try:
                hp2_df = pd.read_csv(loc, sep='\t')
                hp2_df = hp2_df.set_index('# Pathway')
                hp2_df.columns = [s[index_by]]
                hp2_df.index.name = ''
                hp2_df = hp2_df.T
                if norm:
                    count = load_sample(df = s['df'], preproc=s['preproc'], sample = s['fs_name'])
                    comp = general_taxa_one(s)
                    if comp is not None:
                        hp2_df = hp2_df.divide(comp['bacteria'])
                    else:
                        hp2_df = hp2_df.divide(count['reads'])
                dfs.append(hp2_df)
            except:
                print('ERROR: ' + loc)
        else:
            print(loc)
    mp2_combined = pd.concat(dfs, axis=0)
    mp2_combined = mp2_combined.fillna(0)
    
    return mp2_combined

def load_hm2_grouped(samples, index_by='fs_name', norm = False, mapping = 'map_ko_uniref90'):
    dfs = []
    
    db_nucl = 'chocophlan'
    db_protein = 'uniref90'
    
    loc_wc = pipeline + 'datasets/{df}/humann2/{db_nucl}__{db_protein}/{sample}/{preproc}/{sample}__{mapping}.tsv'

    for s in samples:
        loc = loc_wc.format(df = s['df'], 
                            preproc=s['preproc'], 
                            sample = s['fs_name'],
                            db_nucl = db_nucl,
                            db_protein = db_protein,
                            mapping=mapping)
        if os.path.isfile(loc):
            try:
                hp2_df = pd.read_csv(loc, sep='\t')
                hp2_df = hp2_df.set_index('# Gene Family')
                hp2_df.columns = [s[index_by]]
                hp2_df.index.name = ''
                hp2_df = hp2_df.T
                if norm:
                    count = load_sample(df = s['df'], preproc=s['preproc'], sample = s['fs_name'])
                    comp = general_taxa_one(s)
                    if comp is not None:
                        hp2_df = hp2_df.divide(comp['bacteria'])
                    else:
                        hp2_df = hp2_df.divide(count['reads'])
                dfs.append(hp2_df)
            except:
                print('ERROR: ' + loc)
    mp2_combined = pd.concat(dfs, axis=0)
    mp2_combined = mp2_combined.fillna(0)
    
    return mp2_combined

def read_krak_node(df, node_name):
    reads = df.loc[df[5] == node_name][1]
    if len(reads) == 0:
        reads = 0
    else:
        reads = int(reads)
    return reads


def general_taxa_one(s):
    loc_wc = pipeline + 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krak.tsv'
    
    loc = loc_wc.format(df = s['df'], 
                            preproc=s['preproc'], 
                            sample = s['fs_name'],
                            params = 'def')
    
    if os.path.isfile(loc):
        try:
            centr_krak = pd.read_csv(loc, sep='\t', header=None)
            uncl = read_krak_node(centr_krak, 'unclassified')
            vir = read_krak_node(centr_krak, '  Viruses')
            homo = read_krak_node(centr_krak,
                                      '                                                              Homo sapiens')
            bacteria = read_krak_node(centr_krak, '    Bacteria')
            archaea = read_krak_node(centr_krak, '    Archaea')
            other = read_krak_node(centr_krak, 'root') - vir - homo - bacteria - archaea

            total = uncl + vir + bacteria + archaea + homo + other
            composition = {'sample': s['fs_name'],
                               'uncl': uncl,
                               'vir': vir,
                               'bacteria': bacteria,
                               'archaea': archaea,
                               'homo': homo,
                               'other': other,
                               'total': total}
            return composition
        except:
            print('error' + loc)
            return None
    else:
        return None

def get_general_taxa_comp_krak_style(samples):
    loc_wc = pipeline + 'datasets/{df}/taxa/{preproc}/centr__{params}/{sample}/{sample}_krak.tsv'
    comp = []
    
    for s in samples:
        loc = loc_wc.format(df = s['df'], 
                            preproc=s['preproc'], 
                            sample = s['fs_name'],
                            params = 'def')
        if os.path.isfile(loc):
            try:
                centr_krak = pd.read_csv(loc, sep='\t', header=None)
                uncl = read_krak_node(centr_krak, 'unclassified')
                vir = read_krak_node(centr_krak, '  Viruses')
                homo = read_krak_node(centr_krak,
                                      '                                                              Homo sapiens')
                bacteria = read_krak_node(centr_krak, '    Bacteria')
                archaea = read_krak_node(centr_krak, '    Archaea')
                other = read_krak_node(centr_krak, 'root') - vir - homo - bacteria - archaea

                total = uncl + vir + bacteria + archaea + homo + other
                composition = {'sample': s['fs_name'],
                               'uncl': uncl,
                               'vir': vir,
                               'bacteria': bacteria,
                               'archaea': archaea,
                               'homo': homo,
                               'other': other,
                               'total': total}
                comp.append(composition)
            except:
                print('error' + loc)
        else:
            print('NO FILE: ', loc)
    return comp

def load_mag_contigs(meta, source, assembly, centr, binn):
    '''
    Loads info about one bin from MAGs, returns dataframe with contigs coverage info in samples.
    '''
    bin_wc = '/data5/bio/databases/fna/assembly/mh__def/FHM/{ass}/imp__tmtic_def1/conocot_anvio5_def/bin_by_bin/{binn}/{binn}-contigs.names'
    taxa_wc = '/data5/bio/databases/fna/assembly/mh__def/FHM/{ass}/imp__tmtic_def1/conocot_anvio5_def/bin_by_bin/{binn}/{binn}-bin_taxonomy.tab'
    
    samples_for_source = meta.loc[meta['source'] == source]
    samples_for_source = list(samples_for_source['fs_name'])
    T5_stats = bb_stats.get_cov_stats('/data6/bio/TFM/pipeline/datasets', 
                           'FHM', 
                           samples_for_source,
                           'bwa', 
                           'imp__tmtic_def1', 
                           'assembly___mh__def___FHM___{ass}___imp__tmtic_def1'.format(ass = assembly),
                          'final_contigs__1000__no_hum_centr')

    contigs_in_bin = pd.read_csv(bin_wc.format(binn = binn, ass = assembly), header=None)
    contigs_in_bin.columns = ['contig']
    merged = contigs_in_bin.merge(T5_stats, right_on='#ID', left_on='contig')
    merged = merged.drop(['#ID'], axis=1)

    merged['part'] = merged['Length']/merged['Length'].sum()
    for s in samples_for_source:
        merged['avg_on_per__'+s]=merged['Avg_fold__'+s]*merged['Covered_percent__'+s]
#             merged['avg_on_per_on_part__'+s]=merged['avg_on_per__'+s]*merged['part']/meta.loc[meta['fs_name'] == s]['reads'].item()
        merged['avg_on_per_on_part__'+s]=merged['avg_on_per__'+s]*merged['part']/centr.loc[s]['bacteria'].item()
    
    return merged

def load_mags_info(meta, source, assembly, centr):
    '''
    Loads information about MAGs for specific assembly and samples, estimates abundance and returns a dataframe
     with index corresponding to bins and columns corresponding to abundance in samples. Can be transformed to OTU table by applying `df.T`
    '''
    mags = []
    
    # TODO replace with load_mag_contigs function
    bin_wc = '/data5/bio/databases/fna/assembly/mh__def/FHM/{ass}/imp__tmtic_def1/conocot_anvio5_def/bin_by_bin/{binn}'
    taxa_wc = '/data5/bio/databases/fna/assembly/mh__def/FHM/{ass}/imp__tmtic_def1/conocot_anvio5_def/bin_by_bin/{binn}/{binn}-bin_taxonomy.tab'
    bins  = [r.split('/')[-1] for r 
             in glob.glob(bin_wc.format(binn = '*', ass = assembly))]
    bin_wc += '/{binn}-contigs.names'

    samples_for_source = meta.loc[meta['source'] == source]
    samples_for_source = list(samples_for_source['fs_name'])
    T5_stats = bb_stats.get_cov_stats('/data6/bio/TFM/pipeline/datasets', 
                           'FHM', 
                           samples_for_source,
                           'bwa', 
                           'imp__tmtic_def1', 
                           'assembly___mh__def___FHM___{ass}___imp__tmtic_def1'.format(ass = assembly),
                          'final_contigs__1000__no_hum_centr')

    for b in bins:
        contigs_in_bin = pd.read_csv(bin_wc.format(binn = b, ass = assembly), header=None)
        contigs_in_bin.columns = ['contig']
        merged = contigs_in_bin.merge(T5_stats, right_on='#ID', left_on='contig')
        merged = merged.drop(['#ID'], axis=1)

        no_drop = ['contig', 'Length', 'Ref_GC',]

        merged['part'] = merged['Length']/merged['Length'].sum()
        for s in samples_for_source:
            no_drop.append('avg_on_per_on_part__'+s)

            merged['avg_on_per__'+s]=merged['Avg_fold__'+s]*merged['Covered_percent__'+s]
#             merged['avg_on_per_on_part__'+s]=merged['avg_on_per__'+s]*merged['part']/meta.loc[meta['fs_name'] == s]['reads'].item()
            merged['avg_on_per_on_part__'+s]=merged['avg_on_per__'+s]*merged['part']/(centr.loc[s]['bacteria'].item() + centr.loc[s]['uncl'].item())
        cols = list(merged.columns)    
        drop = list(set(cols) - set(no_drop))    

        merged = merged.drop(drop, axis=1)


        taxa = pd.read_csv(taxa_wc.format(binn = b, ass = assembly), header=None, sep='\t')
        taxa = taxa.fillna('Unknown')
        dd = {'Mag': list(taxa[0])[0].split('-')[0] + '__' + list(taxa[1])[0]}
        for s in samples_for_source:
            dd.update({s: merged['avg_on_per_on_part__'+s].sum()})
        mags.append(dd)

    mags = pd.DataFrame(mags)
    mags.index = mags['Mag']
    mags = mags.drop(['Mag'], axis=1)
    return mags