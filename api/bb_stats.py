import pandas as pd
import numpy as np
import os
import yaml
# from Bio import SeqIO
# import matplotlib.mlab as mlab
# import matplotlib.pyplot as plt

config = None
dir_of_this_file = os.path.dirname(os.path.abspath(__file__))
cofig_loc = os.path.join(dir_of_this_file, '../config.yml')
with open(cofig_loc, 'r') as stream:
    try:
        config = yaml.load(stream)
    except yaml.YAMLError as exc:
        print(exc)

ASSNAKE_DB = config.get('assnake_db')
FNA_DB_DIR= config['fna_db_dir']


####------DATA LOADERS------####
def get_cov_stats(prefix, df, samples, tool, preproc, seq_set_path, seq_set_id):
    bb_stats_wc = '{prefix}/{df}/mapped/{tool}__{params}/{path}/{seq_set_id}/{{sample}}/{preproc}/mapped.bb_stats'
    bb_stats_loc = bb_stats_wc.format(
        prefix=prefix,
        df = df,
        preproc=preproc,
        tool = tool,
        params = 'def1',
        path = seq_set_path,
        seq_set_id = seq_set_id
    )

    return load_cov_stats(samples, bb_stats_loc)
    
def load_cov_stats(samples, folder, light = False):
    no_drop = ['#ID', "Length", "Ref_GC", "Avg_fold", 'Norm_fold', "Covered_percent", "Std_Dev"]
    exclude_main = ['Covered_bases',  'Read_GC', 'Plus_reads', 'Minus_reads', 'Median_fold']
    exclude_all = ['Length', 'Ref_GC', 'Covered_bases',  'Read_GC', 'Plus_reads', 'Minus_reads', 'Median_fold']
    if light:
        no_drop =  ['#ID', "Length", "Ref_GC", "Avg_fold", 'Norm_fold', "Covered_percent"]
        exclude_main = ['Covered_bases',  'Read_GC', 'Plus_reads', 'Minus_reads', 'Median_fold', 'Std_Dev']
        exclude_all = ['Length', 'Ref_GC', 'Covered_bases',  'Read_GC', 'Plus_reads', 'Minus_reads', 'Median_fold','Std_Dev']
    dfs = []
    base = folder
    
    for i, s in enumerate(samples):
        fname = base.format(sample=s)
        df = pd.read_csv(fname, sep='\t')
        df['Norm_fold'] = df['Avg_fold']/(df['Plus_reads'] + df['Minus_reads'])
        if i == 0:
            df = df.drop(exclude_main, axis=1)
            df = df[no_drop]
            df = df.add_suffix('__'+s)
            df.rename(columns={'#ID__'+s: '#ID','Length__'+s: 'Length','Ref_GC__'+s: 'Ref_GC' }, inplace=True)
        else:
            df = df.drop(exclude_all, axis=1)
            df = df.add_suffix('__'+s)
            df.rename(columns={'#ID__'+s: '#ID'}, inplace=True)
        dfs.append(df)    
    
    if len(dfs) > 1:
        merged = pd.merge(dfs[0], dfs[1], on ='#ID', how ='inner')
        for i, df in enumerate(dfs):
            if not((i == 0) or (i == 1)):
                merged = pd.merge(merged, dfs[i], on ='#ID', how ='inner')
        return merged
    else: return dfs[0]

def load_coverage(samples, folder, ext='.bb_stats'):
    base = folder
    dfs = []
    for i,s in enumerate(samples):
        cov_file = base+s+ext
        cov_tf = pd.read_csv(cov_file, sep='\t', header = None, names = ['seq', 'start', 'stop', 'cover'])
        dfs.append(cov_tf)
    return dfs

def load_centrifuge(ref):
    cent_f_loc = '/data5/bio/runs-fedorov/TFM/pipeline/data/taxa/assemb/'+ref+'/centr_classification.tsv'
    cent = pd.read_csv(cent_f_loc, sep = '\t')
    
    cent_kr_f_loc = '/data5/bio/runs-fedorov/TFM/pipeline/data/taxa/assemb/'+ref+'/centr_krak.tsv'
    #(U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, (S)pecies).
    cent_krak = pd.read_csv(cent_kr_f_loc, sep = '\t', header = None)
    cent_krak.columns = ['percent', 'num_of_reads', 'num_no_futher_class', 'code', 'NCBI_id', 'name']
    
    cent_cl = '/data5/bio/runs-fedorov/TFM/pipeline/data/taxa/assemb/'+ref+'/centr_report.tsv'
    cent_cl = pd.read_csv(cent_cl, sep = '\t')
    return cent, cent_krak, cent_cl
  
############----PLOTTING----###########
def to_nucl_res(data):
    nuc = np.empty(int(data.iloc[-1].stop))
    for (i, s,t,c) in data.itertuples():
        nuc[s:t] = c
    return pd.DataFrame(dict(start = range(0,len(nuc)), cover = nuc))

def roller(data, window):
    d = []
    for i in range(0, len(data), window):       
        d.append([i, data[i: i + window].cover.mean()])
    return pd.DataFrame(d, columns=['start', 'cover'])

def prepare(data, org, norm, roll = False, remove_outliers = False, agg_win = 100):
    # agg_win = something like seq_len/4000
    agg_win = 2
    roll_win = 10
    subset = data[data['seq'] == org][['start', 'stop', 'cover']]
    subset.cover = subset.cover/float(norm)
    d = to_nucl_res(subset)
    if remove_outliers:
        d = remove_outliers(d, 'cover')
    if roll:
        d = roller(d, agg_win)
        d = d.rolling(window = roll_win, center = True, min_periods = 1, on = 'start').mean()
    return d

def remove_outliers(df_in, col_name):
    q1 = df_in[col_name].quantile(0.25)
    q3 = df_in[col_name].quantile(0.75)
    iqr = q3-q1 #Interquartile range
    fence_low  = q1-1.5*iqr
    fence_high = q3+1.5*iqr
    #df_in[col_name] = df_in[col_name] < fence_high
    df_in.loc[df_in.cover > fence_high, 'cover'] = fence_high
    return df_in

def plot_coverage(seq, samples, cov_dfs_list, cov_info_all, roll = False, remove_outliers = False):
    proc = []
    seq_cov_info = cov_info_all.loc[cov_info_all['#ID'] == seq]
    length = int(seq_cov_info['Length'])
    for i, s in enumerate(samples):
        curr_cov_df = cov_dfs_list[i]
        norm = float(seq_cov_info['Avg_fold__'+s]/seq_cov_info['Norm_fold__'+s])
        proc.append(prepare(cov_dfs_list[i], seq, norm, roll, remove_outliers, int(length/4000)))
        
    fig = plt.figure(1, figsize = (12,6))
    plt.xlabel('Nucleotide')
    plt.ylabel('Coverage')
    plt.title(seq)
    for i, s in enumerate(samples):
        plt.plot(proc[i].start, proc[i].cover, label=s)
    plt.legend()
    plt.show()

#########################

# Used to print sequnce from fasta by id
def print_seq(seq_id, seq_dict):
    print ('>'+seq_id)
    print (seq_dict[seq_id].seq)

#used to plot portrait
def plot_gc_cov_portrait(gc, cov, ax = None, rgba_colors = [[1,1,1,1]], size=10, title='GC-cov-taxa Portrait'):
    fig, ax = plt.subplots(figsize=(size,size))
    ax.scatter(x = gc, 
    y = cov, 
    # color = rgba_colors
    )
    ax.set_title(title)
    ax.set_xlabel('GC content')
    ax.set_ylabel('log cov')

    plt.show() 

#many portraits
def plot_gc_cov_portrait_mult(df, samples, ax, norm = False, size = 4, title='cov-gc poratrait'):
    for s in samples:
        yy = np.log(df['Avg_fold__'+s]+1)
        if norm:
            yy = np.log(df['Norm_fold__'+s]+1)
        ax.scatter(x = df['Ref_GC'], y = yy, label = s)
    ax.set_title(title)
    ax.set_xlabel('GC content')
    ax.set_ylabel('log cov')
    ax.legend()
    
def plot_portrait_with_diff(omg, sample, s_vs_s, ax, size = 4):
    alphas = []
    colors = []
    rgba_colors = np.zeros((len(omg),4))
    same = 0
    for i, per in enumerate(omg[s_vs_s]):
        aa = 0
        if per > 98:
            rgba_colors[i,0] = 1.0
            rgba_colors[i, 3] = 1.0
            same+=1
        elif ((per < 98) & (per > 80)) :
            rgba_colors[i,0] = 0.3
            rgba_colors[i,1] = 1.0
            rgba_colors[i, 3] = 0.75
        elif ((per < 80) & (per > 65)) :
            rgba_colors[i,2] = 1.0
            rgba_colors[i, 3] = 0.47
        else :
            rgba_colors[i,0] = 0
            rgba_colors[i,1] = 0
            rgba_colors[i,2] = 0
            rgba_colors[i, 3] = 0.15
        alphas.append(aa)
    plot_gc_cov_portrait(omg.Ref_GC, omg[sample], ax, rgba_colors, 
                         size = size, title = (s_vs_s[0:-6] + '; N_same: '+str(same)))
    
def get_df_from_query(data, q):
    # Because it doesn't accept these symbols in query, remove
    delete_symb = str.maketrans(dict.fromkeys("-_"))
    q = q.translate(delete_symb)
    
    df = data
    cols_before = data.columns
    cols = cols_before.map(lambda x: x.replace('_', '').replace('-', ''))
    df.columns = cols
    df = df.query(q)
    df.columns = cols_before
    data.columns = cols_before
    
    return df

def basic_info(data, samples, queries, titles, print_stats = False):
    fig, axes = plt.subplots(nrows=1, ncols=len(queries), figsize=(15, 6))
    
    for i, q in enumerate(queries):
        df = get_df_from_query(data,q)
        plot_gc_cov_portrait_mult(df, samples, axes[i], norm=True, title=titles[i]+'; N = ' + str(len(df)))
    
    fig.tight_layout()
    plt.show()
    if print_stats:
        for s in samples:
            print('Mean cov % in ' + s +': ' + str(data['Covered_percent__'+s].mean()))
            print('Std Dev cov % in ' + s +': ' + str(data['Covered_percent__'+s].std()))
            print('Median cov % in ' + s +': ' + str(data['Covered_percent__'+s].median()) +'\n')