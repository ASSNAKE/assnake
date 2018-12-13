#/home/fedorov/miniconda3/bin python 

import os, sys, getopt

#Assumption - It reads ouput of multisample SNP calling from GATK. The output is parsed using VCFTools to include just two samples for comparision. The samples two lines of input looks like as giebn below
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  D1      RpT1m6
#Bacteroides_vulgatus_PC510      784     .       A       G       811.11  .       .       GT:AD:DP:GQ:PL  0:47,0:47:99:0,2004     0:79,0:79:99:0,3361

# read input file as command line option

infile=''
outfile=''

# defines a function called checkwindow which input two arguments like col9 and col10, which are two samples of a multiple vcf file and return two values (if they have same SNP pattern and how many SNPs are found in total. This ignore SNPs if they any sample has allele freq >25%.

def checkwindow(col9,col10):
   
    # Allele frequency cutoff for 20%
    ALLEL_FREQ1 = 0.20
    ALLEL_FREQ2 = 5.0
    LOW_COV = 4

    #grab GT
    data_col9 =  col9.split(':')
    data_col10 =  col10.split(':')
    
    #grab AD and their count
    list2_col9 = data_col9[1].split(',')
    list2_col10 = data_col10[1].split(',')


    sum_list2_col9 = 0
    sum_list2_col10 = 0
    if data_col9[1] == '.':
        sum_list2_col9 = 0
    else:
        sum_list2_col9 = sum(map(int,list2_col9)) 

        if data_col10[1] == '.':
            sum_list2_col10 = 0
        else:
            sum_list2_col10 = sum(map(int,list2_col10))

    min_depth = 0
    #calculate minimum depth of two above
    if sum_list2_col9 <= sum_list2_col10:
        min_depth = sum_list2_col9
        
    else:
        min_depth = sum_list2_col10
    
    # Check and mark ignore  multiallele sites
    if data_col9[0] == '.' or data_col10[0] == '.' or min_depth <= LOW_COV:
        #print "MULTI ALLELE ",data_col9[0] ,data_col10[0]
        return ("LOW_COV:DEP-" + str(min_depth))

    # Check and mark ignore  multiallele sites
    elif int(data_col9[0]) > 1 or int(data_col10[0]) > 1:
        #print "MULTI ALLELE ",data_col9[0] ,data_col10[0]
        return "SKIP_MULTIALLELE:DEP-" + str(min_depth)

    # Checking if any SNP present and if they have allel frequency > 20
    elif int(data_col9[0]) == 1 or int(data_col10[0]) == 1:
        #print data_col9[0] ,data_col10[0]
        # to calculate  allele frequency
                #print list2_col9[0],"``",list2_col9[1]
                #print list2_col10[0],"``",list2_col10[1]

        # accomodating for division by zero error
 
        try:
            div_col9 = float(list2_col9[0]) /float(list2_col9[1])
        except ZeroDivisionError:
            div_col9 = 0
        
        try:
            div_col10 = float(list2_col10[0]) /float(list2_col10[1])
        except ZeroDivisionError:
            div_col10 = 0

        if (div_col9 <= ALLEL_FREQ1 or div_col9 >= ALLEL_FREQ2) and (div_col10 <= ALLEL_FREQ1 or div_col10 >= ALLEL_FREQ2):
            if int(data_col9[0]) == 1 and int(data_col10[0]) == 1:
                return "SKIP_SAME_SNP:DEP-" + str(min_depth)
            elif int(data_col9[0]) == 0 and int(data_col10[0]) == 1:
                return "DIFFERENT_SNP_01:DEP-" + str(min_depth)
            else:
                return "DIFFERENT_SNP_10:DEP-" + str(min_depth)
        else:
            return "SKIP_HighAlleleFreq:DEP-" + str(min_depth)

    elif int(data_col9[0]) == 0 and int(data_col10[0]) == 0:
        return "SKIP_NO_SNP:DEP-" + str(min_depth)

    else:
        return "ERROR:UNEXPECTED OPTION:DEP-" + str(min_depth)



# defines arguments and how to process them
myopts, args = getopt.getopt(sys.argv[1:],"i:o:")

for o, a in myopts:
    if o == '-i':
        infile=a
    elif o == '-o':
        outfile=a
    else:
        print("Usage: %s -i input -o output" % sys.argv[0])


#print ("Input file : %s and output file: %s" % (infile,outfile) )


#Open input file
INFILE = ''
OUTFILE = ''
if os.path.isfile(infile):
    try:
        INFILE  = open(infile, "r")
        print ("\nINFO: Reading input file " + infile)
    except IOError:
        print ('ERROR: Cannot open file ', + infile )
        sys.exit()
else:
    print ("ERROR: Input file " + infile + " not found")
    sys.exit()

try:
    OUTFILE  = open(outfile, "w")
    print ("INFO: Opened  output file " + outfile + " for writing")
except IOError:
    print ('ERROR: Cannot open file ', + outfile)
    sys.exit()



count_TOTAL = 0
count_OTHER = 0
count_SKIP_MULTIALLELE = 0
count_SKIP_SAME_SNP = 0
count_DIFFERENT_SNP_01 = 0
count_DIFFERENT_SNP_10 = 0
count_SKIP_HighAlleleFreq = 0
count_SKIP_NO_SNP = 0
count_LOW_COV = 0

print ("INFO: Analysing data, please wait...")

# with open(infile) as f:
#     for line in f:
#         if not line.startswith('##contig'): 
#             print(line)
# print('Thats all, folks!')

for line in INFILE:
#     print (line)
    # Skipping the lines starting with '#'
    if not line.startswith('#'):
        #print line
        line = line.rstrip('\r\n|\n')  # to remove endlines if any
        info = line.split()
        try:
            col9 = info[9]
        except IndexError:
            print ('ERROR: SNP data for first sample not found. Terminating program ...')
            OUTFILE.close()
            os.remove(outfile)
            sys.exit()

        try:
            col10 = info[10]
        except IndexError:
            print ('ERROR: SNP data for second sample not found. Terminating program ...')
            OUTFILE.close()
            os.remove(outfile)
            sys.exit()
        
        result = checkwindow(col9, col10)
        OUTFILE.write(line + "\t" + result + "\n")

        #count different types of SNPs comparision
        result = result.split(':')[0]

        if result == 'SKIP_NO_SNP':
            count_SKIP_NO_SNP += 1
        elif result == 'SKIP_HighAlleleFreq':
            count_SKIP_HighAlleleFreq += 1
        elif result == 'DIFFERENT_SNP_01':
            count_DIFFERENT_SNP_01 += 1
        elif result == 'DIFFERENT_SNP_10':
            count_DIFFERENT_SNP_10 += 1
        elif result == 'SKIP_SAME_SNP':
            count_SKIP_SAME_SNP += 1
        elif result == 'SKIP_MULTIALLELE':
            count_SKIP_MULTIALLELE += 1
        elif result == 'LOW_COV':
            count_LOW_COV += 1
        else:
            count_OTHER += 1

        count_TOTAL += 1
        
        #print line,"\t",result 
        #print line

    else:
        OUTFILE.write(line)

print ("INFO: Output written in file named " + outfile)
print ("\nINFO: -------------------STATS----------------------")
print ("INFO: Total-sites\t", count_TOTAL)
print ("INFO: Ignored low coverage sites\t", count_LOW_COV)
print ("INFO: Ignored Multiallelic sites\t",count_SKIP_MULTIALLELE)
print ("INFO: Ignored HighAllele frequency SNPs", count_SKIP_HighAlleleFreq)
print ("INFO: No SNPs observed\t", count_SKIP_NO_SNP)
print ("INFO: Same SNPs observed\t", count_SKIP_SAME_SNP)
print ("INFO: Different SNPs observed as 0 1\t", count_DIFFERENT_SNP_01)
print ("INFO: Different SNPs observed as 1 0\t", count_DIFFERENT_SNP_10)
print ("INFO: -------------------STATS----------------------\n")
print ("INFO: Program finished sucessfully ...\n")
OUTFILE.close()
INFILE.close()