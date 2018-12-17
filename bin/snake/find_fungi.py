import os
rule fungi_kraken:
    input:
        FindFungi1 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi1/{sample}/report.tsv',
        FindFungi2 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi2/{sample}/report.tsv',
        FindFungi3 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi3/{sample}/report.tsv',
        FindFungi4 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi4/{sample}/report.tsv',
        FindFungi5 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi5/{sample}/report.tsv',
        FindFungi6 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi6/{sample}/report.tsv',
        FindFungi7 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi7/{sample}/report.tsv',
        FindFungi8 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi8/{sample}/report.tsv',
        FindFungi9 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi9/{sample}/report.tsv',
        FindFungi10 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi10/{sample}/report.tsv',
        FindFungi11 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi11/{sample}/report.tsv',
        FindFungi12 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi12/{sample}/report.tsv',
        FindFungi13 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi13/{sample}/report.tsv',
        FindFungi14 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi14/{sample}/report.tsv',
        FindFungi15 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi15/{sample}/report.tsv',
        FindFungi16 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi16/{sample}/report.tsv',
        FindFungi17 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi17/{sample}/report.tsv',
        FindFungi18 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi18/{sample}/report.tsv',
        FindFungi19 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi19/{sample}/report.tsv',
        FindFungi20 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi20/{sample}/report.tsv',
        FindFungi21 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi21/{sample}/report.tsv',
        FindFungi22 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi22/{sample}/report.tsv',
        FindFungi23 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi23/{sample}/report.tsv',
        FindFungi24 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi24/{sample}/report.tsv',
        FindFungi25 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi25/{sample}/report.tsv',
        FindFungi26 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi26/{sample}/report.tsv',
        FindFungi27 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi27/{sample}/report.tsv',
        FindFungi28 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi28/{sample}/report.tsv',
        FindFungi29 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi29/{sample}/report.tsv',
        FindFungi30 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi30/{sample}/report.tsv',
        FindFungi31 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi31/{sample}/report.tsv',
        FindFungi32 = 'datasets/{df}/taxa/{preproc}/kraken-v1.1-ff/FindFungi32/{sample}/report.tsv',
    output: all_class_sorted = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/all_classified_sorted.tsv'
    params:
        all_class = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/all_classified.tsv',
        all_class_sorted = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/all_classified_sorted.tsv',
        read_names = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/unique_reads.txt',
        reads_fasta = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/unique_reads.fa',
        reads_fasta_ref = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/unique_reads_reformatted.fa',
    threads: 6
    run:
        all_sorted_str = ''

        for report in input:
            print('filtering ' + report.split('/')[-3])
            shell('grep ^C {report} >> {params.all_class}')

            report_classified = os.path.dirname(report) + '/report_classified.tsv'
            report_classified_sorted = os.path.dirname(report) + '/report_classified_sorted.tsv'
            shell('grep ^C {report} > {report_classified}')
            # Sort individual Kraken output files
            shell('sort -o {report_classified_sorted} -k2,2 {report_classified}')
            all_sorted_str += report_classified_sorted + ' '

        # Merge and sort all Kraken output files
        shell('sort -o {output.all_class_sorted} -m -k2,2 {all_sorted_str}')

        shell("cat {output.all_class_sorted} | awk '{{print $2}}' | sort | uniq > {params.read_names}")

        unique_reads = []
        with open(params.read_names, 'r') as unique:
            for line in unique.readlines():
                unique_reads.append(line[0:-1])
        print(len(unique_reads))
        saved_read_names = set([])
        reads_num = 0

        # gather reads in one fasta
        for report in input:
            unqiue_reads_str = ''
            reads = os.path.dirname(report) + '/classified.tsv'
            print(reads)
            with open(reads, 'r') as classified_reads:
                cnt = 1
                line = classified_reads.readline()
                while line:
                    if cnt%2 == 1:
                        read = line[0:-1].split(' ')[0]
                        if read not in saved_read_names:
                            saved_read_names.add(read)
                            unqiue_reads_str += read + '\n'
                            unqiue_reads_str += classified_reads.readline()
                            cnt += 1
                        else:
                            cnt += 1
                            classified_reads.readline()
                        line = classified_reads.readline()
                        cnt += 1
            with open(params.reads_fasta, 'a') as reads_fasta:
                reads_fasta.write(unqiue_reads_str)
        shell('''awk 'NR==1{{printf $0"\\t";next}}{{printf /^>/ ? "\\n"$0"\\t" : $0}}' {params.reads_fasta} > {params.reads_fasta_ref}''')
        print(len(saved_read_names))
        # shell ('touch {output}')

consensus_script = config['FindFungi']['v0.23.3']['consensus']
taxids = config['FindFungi']['v0.23.3']['taxids']
rule fungi_step2:
    input: all_class_sorted = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/all_classified_sorted.tsv'
    output: consensus = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/consensus.tsv',
            predictions =    'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/predictions.tsv',
    params: wd =    'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/processing/',
    conda: "/data6/bio/TFM/pipeline/envs/FindFungi.yml"

    shell: ("""python2.7 {consensus_script} {input.all_class_sorted} {output.consensus}; \
        awk '{{print $3}}' {output.consensus} | sort -n | uniq -c | sort -k1,1nr > {output.predictions}; \
        echo "Generated predictions" \n
        mkdir -p {params.wd} \n
        while read p; do \n
            ReadNumber=$(echo $p | awk -F ' ' '{{print $1}}') \n
            Taxid=$(echo $p | awk -F ' ' '{{print $2}}') \n
            if [ $ReadNumber -ge 10 ] \n
            then \n
                if grep -Fxq $Taxid {taxids} \n
                    then #If the predicted taxid is one of the 949 fungal species \n
                        echo "BLASTing $Taxid" \n
                        #Gather read names for each taxid \n
                        mkdir -p {params.wd}$Taxid
                        awk -v taxid="$Taxid" '$3 == taxid {{print $2}}' {output.consensus} > {params.wd}$Taxid/read_names.txt \n
                    else \n
                        echo "Not BLASTing $Taxid, wrong taxonomy level" \n
                fi \n
            else \n
                echo "Not BLASTing $Taxid, too few reads" \n
            fi \n
        done < {output.predictions}""")

rule ff_single_blast:
    input:
        read_names = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/processing/{taxid}/read_names.txt',
        reads_fasta_ref = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/unique_reads_reformatted.fa',
    output:
        # blast = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/processing/{taxid}/blast.tsv',
        fasta = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/processing/{taxid}/reads.fa'
    run:
        print('Gathering reads for taxid ' + wildcards.taxid)
        shell('''awk -v reads="{input.read_names}" -F "\\t" 'BEGIN{{while((getline k < reads)>0)i[k]=1}}{{gsub("^>","",$0); if(i[$1]){{print ">"$1"\\n"$2}}}}' {input.reads_fasta_ref} > {output.fasta}''')

def get_taxids(wildcards):
    dir = 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/processing/*'
    dir = dir.format(df = wildcards.df, preproc=wildcards.preproc, sample=wildcards.sample)
    taxids = [t+'/reads.fa' for t in glob(dir)]
    return taxids

rule find_fungi_blastn:
    input: get_taxids
    output: 'datasets/{df}/taxa/{preproc}/FindFungi/{sample}/blast.done',
    run:
        for f in input:
            print(f)
