# from Bio import SeqIO
# import pysam

def load_coverage_info_for_markers_in_sequences(fasta_loc, bampath_wc, working_dir,
                                                only_variable = True, print_progress = True):
    fasta_sequences = list(SeqIO.parse(open(fasta_loc),'fasta'))

    vrb = []
    cov_dfs = []
    bampath = bampath_wc
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
        
    with pysam.AlignmentFile(bampath, 'rb') as bamfile:
        for sequence in fasta_sequences:
            
            if print_progress: print(sequence.id)

            start = 0
            end = len(sequence.seq)

            counts = bamfile.count_coverage(
                sequence.id, 
                start=start, 
                end=end, 
                quality_threshold=10)
            seq_wd = os.path.join(working_dir, sequence.id, 'sequence'+'.cov')
            os.mkdir(os.path.join(working_dir, sequence.id))
            lines = ['ref_allele\tref_pos\tdepth\tA\tT\tG\tC\n']
            for i in range(start, end):
                depth = sum([counts[_][i] for _ in range(4)])

                if depth > 0:
                    ref_pos = i+1
                    ref_allele = sequence.seq[i]
                    a = counts[0][i]
                    c = counts[1][i]
                    g = counts[2][i]
                    t = counts[3][i]
                        
                    lines.append('{ref_allele}\t{ref_pos}\t{depth}\t{A}\t{T}\t{G}\t{C}\n'.format(ref_pos= ref_pos, 
                        ref_allele= ref_allele, 
                        depth= depth, 
                        A= a, T= t, G=g, C= c,
                    ))
                
            with open(seq_wd, 'x') as cov_file:
                cov_file.writelines(lines)         
    return cov_dfs

rule produce_cov:
    input: 
        bam = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_kumar.bam',
        fa = get_reference_fasta
    output:
        done = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_kumar__msnv/msnv.done'
    params: wd = '{prefix}/{df}/mapped/bwa__{params}/{path}/{seq_set_id}/{sample}/{preproc}/mapped_kumar__msnv/'
    run:
        load_coverage_info_for_markers_in_sequences(input.fa, input.bam, params.wd)
        shell('touch {output.done}')