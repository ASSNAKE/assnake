
def before(value, a):
    # Find first part and return slice before it.
    pos_a = value.find(a)
    if pos_a == -1: return ""
    return value[0:pos_a]

def get_samples_for_bowtie2(wildcards):
    folder = "datasets/"+wildcards.df+"/reads/"+wildcards.preproc+'/'
    sample_id = wildcards.id_sample
    files = []
    result = ''
    res = []
    
    for fname in os.listdir(folder):
        fl = before(fname, ".fastq")
        last_3 = (fl[len(fl)-3:len(fl)])
        if last_3 == '_R1' or last_3 == '_R2':
            fl = fl[0:len(fl)-3]
        if sample_id == fl:
            files.append(fname)
    i = 0
    if len(files) < 1:
        raise Exception("NO FILES!")
    elif len(files) > 1:
        for f in files:
            if '_S.' not in f:
                result += f
                res.append(folder+f)
                if i == 0:
                    result += ' '
                    i += 1
    elif len(files) == 1:
        res.append(folder+str(files[0]))
        result += files[0]
    return res


#Creates sequence index for bowtie2
rule create_seq_set_index_bowtie2:
    input:
        ref = "data/ref/{type}/{seq_set_id}.fa"
    output:
        _1bt2 = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.1.bt2',
        _2bt2 = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.2.bt2',
        _3bt2 = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.3.bt2',
        _4bt2 = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.4.bt2',
        ref1bt2 = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.rev.1.bt2',
        ref2bt2 = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.rev.2.bt2'
    params:
        prefix = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}'
    log: 'log/bowtie2/index/{type}/{seq_set_id}.log'
    run:
        shell('echo -e "INFO: Creating Bowtie2 index for {input.ref}\n"')
        shell('''bowtie2-build {input.ref} {params.prefix} > {log} 2>&1 ''')
        shell('echo -e "INFO: Finished creating BWA index for {input.ref}\n"')

        
rule bowtie2:
    input:
        r1 = 'datasets/{df}/reads/{preproc}/{sample_id}_R1.fastq.gz',
        r2 = 'datasets/{df}/reads/{preproc}/{sample_id}_R2.fastq.gz',
        s = 'datasets/{df}/reads/{preproc}/{sample_id}_S.fastq.gz',
        ref_index = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}.1.bt2'
    output:
        sam = 'datasets/{df}/mapped/{preproc}/bowtie2/{type}/{seq_set_id}/mapped/{sample_id}.sam'
    params: 
        index = 'data/ref/index/bowtie2/{type}/{seq_set_id}/{seq_set_id}'
    log: 'log/bowtie2/{df}/reads/{preproc}/{type}/{seq_set_id}/{sample_id}.log'
    threads: 2
    run:
        shell('bowtie2 -p {threads} --mm -x {params.index} --no-unal -1 {input.r1} -2 {input.r2} -U {input.s} -k 1 -S {output.sam} >{log} 2>&1')
        
    