import pandas as pd

def create_depth_table_for_metabat2(samples):
    print('in func')
    print(list(samples))
    print(type(samples))

    samples = list(samples)
    col_list = []
    # read first sample
    total_depth = pd.read_csv(samples[0], sep='\t')
    s_name = samples[0].split('/')[-3]
    rename_dict = {'mapped.bam':'{s}.bam'.format(s=s_name), 'mapped.bam-var':'{s}.bam-var'.format(s = s_name)}
    col_list.append('{s}.bam'.format(s=s_name))
    total_depth = total_depth.rename(rename_dict, axis=1)
    total_depth = total_depth.drop('totalAvgDepth', axis=1)
    
    for sample in samples[1:]:
        depth = pd.read_csv(sample, sep='\t')
        s_name = sample.split('/')[-3]
        col_list.append('{s}.bam'.format(s=s_name))
        rename_dict = {'mapped.bam':'{s}.bam'.format(s=s_name), 'mapped.bam-var':'{s}.bam-var'.format(s = s_name)}
        depth = depth.rename(rename_dict, axis=1)
        depth = depth.drop(['totalAvgDepth', 'contigLen'], axis=1)
        total_depth = total_depth.merge(depth, on = 'contigName')
    
    totalAvgDepth = total_depth[col_list].sum(axis=1)
    total_depth.insert(2,'totalAvgDepth',totalAvgDepth)
    return total_depth

def get_samples_for_metabat(wildcards):
    bam_wc = '{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/mapped.jgi_depth'
    ss = wildcards.samples

    table_wc = '{prefix}/{df}/assembly/mh__{params}/{samples}/sample_set.tsv'
    r_wc_str = '{prefix}/{df}/reads/{preproc}/{sample}/{sample}_{strand}.fastq.gz'
    
    table = pd.read_csv(table_wc.format(prefix = wildcards.prefix,
                                        df = wildcards.df,
                                        params = wildcards.params,
                                        samples = wildcards.samples),
                        sep = '\t')

    list_of_sample_profiles = []
    for s in table.to_dict(orient='records'):     
        list_of_sample_profiles.append(bam_wc.format(prefix=wildcards.prefix,
                                    bwa_params = wildcards.bwa_params,
                                    params = wildcards.params,
                                    mod =  wildcards.mod, 
                                    samples = wildcards.samples,
                                    df=s['df'], 
                                    preproc=s['preproc'],
                                    sample=s['fs_name']))

    print(list_of_sample_profiles)
    return list_of_sample_profiles

rule jgi_sum_depth:
    input:'{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/mapped.bam'
    output:'{prefix}/{df}/mapped/bwa__{bwa_params}/assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}/{sample}/{preproc}/mapped.jgi_depth'
    conda: 'metabat2_env.yaml'
    shell: ('jgi_summarize_bam_contig_depths --outputDepth {output} {input}')

rule depth_file:
    input: get_samples_for_metabat
    output: '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/depth.tsv'
    run: 
        depth = create_depth_table_for_metabat2(input)
        print(output)
        depth.to_csv(str(output), sep='\t', index=False)

rule metabat2_new:
    input:
        fa    = os.path.join(fna_db_dir, 'assembly/mh__{params}/{df}/{samples}/final_contigs__{mod}.fa'),
        depth = '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/depth.tsv'
    output:
        done  = '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/metabat2.done'
    params:
        wd    = '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/',
        bin_d = '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/bins/bin',
        bin_b = '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/bins/',
    log:        '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/log.txt'
    benchmark:  '{prefix}/{df}/metabat2/bwa__{bwa_params}/mh__{params}/{samples}/{mod}/benchmark.txt'
    threads: 20
    conda: 'metabat2_env.yaml'
    shell: ('''cd {params.wd}; \n
            (metabat2 -t {threads} -i {input.fa} -a {input.depth} -o {params.bin_d}) >{log} 2>&1; \n
            cd {params.bin_b}; \n
            for f in *.fa; do i="${{f%.fa}}"; mv $f ${{i//./_}}.fa; done; \n
            touch {output.done}''')

def external_binning():
    bins_wc = '/data10/bio/metagenomics/FHM/metabat2/bwa__def/mh__def/D2/1000/bins/{binn}.fa'
    bins = [b.split('/')[-1][0:-3] for b in glob.glob(bins_wc.format(binn = '*'))]

    bb = []
    bins_info = []
    for i, b in enumerate(bins):
        with open(bins_wc.format(binn = b), 'r') as contigs:
            for line in contigs:
                if line[0] == '>':
                    bb.append({ 'contig': line[1:-1], 'bin':b.replace('.', '_')})
                    
    external_binning_of_contigs = '/data10/bio/metagenomics/FHM/metabat2/bwa__def/mh__def/D2/1000/external_binning_of_contigs.tsv'
    external_binning = pd.DataFrame(bb, columns=['contig', 'bin'])
    external_binning.to_csv(external_binning_of_contigs, index = None, header=False, sep='\t')