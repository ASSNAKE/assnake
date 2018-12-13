bowtie2 -U kazan5500_2014_07_10_4_L01_GorK_F3.fastq -x /data3/bio/human_genome/bowtie2_hg19/hg19 --no-unal kazan5500_2014_07_10_4_L01_GorK_F3_nohuman.fastq -p 20 1> /dev/null 2> bowtie2.log" | qsub -cwd -pe make 20

