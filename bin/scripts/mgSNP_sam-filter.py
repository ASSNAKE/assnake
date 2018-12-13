#Note from KUMAR
#Note I noticed that some program modify the column number of MD data. 
#please check if MD column used here is 11 (starting from 0). Otherwise change line 121

#/home/fedorov/miniconda3/bin python 
import os, sys, getopt, re

MIN_ALIGN_LENGTH = 80
MIN_PERCENT_MATCH = 92

# Here tab5 is the CIGAR string
def getCIGARstats(tab5):
    readLength = 0
    clipLength = 0
    indelLength = 0
    #print tab9,
    #seqlength = len(tab9)
    #print seqlength,
    #print tab5,
    tab5 = tab5.strip()
    match_all = re.findall( r"([0-9]*[MIDNSHP=X]*)", tab5, re.I)
    # ignoring tag except H,S,M
    for item in match_all:
        #print item,
        if re.match('.*H|.*S',item):
            tmp_item = item.strip('HS')
            clipLength += int(tmp_item)
            readLength += int(tmp_item)
            #print "-H/S-",
        elif re.match('.*M|.*=|.*X',item):
            tmp_item = item.strip('MX=')
            readLength += int(tmp_item)
            #print "-M/=/X-",
        elif re.match('.*I|.*D',item):
            tmp_item = item.strip('ID')
            indelLength += int(tmp_item)
            #print "-I/D-",
        else:
            pass
    #print "*", readLength, "*",clipLength,"*", indelLength,"*"
    return readLength, clipLength, indelLength

def getMDstats(tab11):
    #MD = "0G24A1T22^CTA72T53^A87A76"
    #print "MD:" + MD
    MD_tag = tab11.split(':')
    MD = MD_tag[2]
    matchLength = 0
    mismatchLength = 0
    deletionLength = 0

    match_all = re.findall( "(\^[A-Z]*)", MD, re.I)
    for item in match_all:
        if item != '':
            #print item,
            tmp_item = item.strip('^')
            deletionLength += len(tmp_item)

    match_all = re.findall( "([A-Z]*)", MD, re.I)
    for item in match_all:
        if item != '':
            #print item,
            mismatchLength += len(item)
    #print

    match_all = re.findall( "([0-9]*)", MD, re.I)
    for item in match_all:
        if item != '':
            #print item,
            matchLength += int(item)
        #print

    mismatchLength = mismatchLength - deletionLength
    #print "Match-" + str(matchLength) + "\tMisMatch-" + str(mismatchLength) + "\tDeletion-" + str(deletionLength)
    return matchLength, mismatchLength, deletionLength

def filter_sam(INFILE, OUTFILE):
    i = 0
    pr = True
    for line in INFILE:
        #print line
        # Skipping the lines starting with '#'
        if not line.startswith('@'):
            #print(line)
            line = line.rstrip('\r\n|\n')  # to remove endlines if any
            tab = line.split()
            tab2 = tab[2]
            tab5 = tab[5]
            tab9 = tab[9]
            tab12 = tab[12]
            if tab2 == '*' or tab5 == '*':
                #OUTFILE.write(line + "\n")
                pass
            else:
                if pr:
                    print (tab12)
                    pr = False
                readLength, clipLength, indelLength = getCIGARstats(tab5)
                alignLength = readLength - clipLength
                matchLength, mismatchLength, deletionLength = getMDstats(tab12)
                percentMatch = 100 * (float(matchLength) / float(matchLength + mismatchLength))
                #print "*", readLength, "*",clipLength,"*", indelLength,"*"
                if alignLength > MIN_ALIGN_LENGTH and percentMatch >= MIN_PERCENT_MATCH:
                    #print line
                    OUTFILE.write(line + "\n")

                    #print(tab5,tab12 , matchLength, mismatchLength, deletionLength, alignLength , percentMatch)

        else:
            OUTFILE.write(line)

def main():
    #default vars
    MIN_ALIGN_LENGTH = 80
    MIN_PERCENT_MATCH = 92
    infile = ''
    outfile = ''
    
    # defines arguments and how to process them
    myopts, args = getopt.getopt(sys.argv[1:],"i:o:m:l:")
    for o, a in myopts:
        if o == '-i':
            infile=a
        elif o == '-o':
            outfile=a
        elif o == '-m':
            MIN_PERCENT_MATCH = int(a)
        elif o == '-l':
            MIN_ALIGN_LENGTH = int(a)
        else:
            print("Usage: %s -i input -o output -m percent_match -l min_al_len" % sys.argv[0])
    
    #Open input file
    if os.path.isfile(infile):
        try:
            INFILE  = open(infile, "r")
            print("\nINFO: Reading input file " + infile)
        except IOError:
            print('ERROR: Cannot open file ' + infile)
            sys.exit()
    else:
        print("ERROR: Input file " + infile + " not found")
        sys.exit()

    #Open output file
    try:
        OUTFILE  = open(outfile, "w")
        print("INFO: Opened  output file " + outfile + " for writing")
    except IOError:
        print('ERROR: Cannot open file ' + outfile)
        sys.exit()

    filter_sam(INFILE, OUTFILE)
    
    print ("INFO: Program completed. Output file  " + outfile + " generated.\n")

#print ("Input file : %s and output file: %s" % (infile,outfile) )

if __name__ == "__main__":
    main()




