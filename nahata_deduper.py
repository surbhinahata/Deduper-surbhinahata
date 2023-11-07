#!/usr/bin/env python
# author: Surbhi Nahata
# Date: 11/2/2023

###########
# imports #
###########
import argparse #argparse to make the code general for any file
import re #regex expression for adjusted position

############################
# arguments and functions #
############################

def get_args():
    parser = argparse.ArgumentParser(description="a function to remove PCR duplicates from sequencing data") 
    parser.add_argument("-f", "--file", help = "designates absolute file path to sorted sam file", required = True)
    parser.add_argument("-o", "--outfile", help = "designates absolute file path to sorted sam file", required = True)
    parser.add_argument("-u", "--umi", help = "designates file containing the list of UMIs", required = True)
    args = parser.parse_args()
    return args
args = get_args()

infile = args.file #input file
outfile = args.outfile #output file with 
umifile = args.umi

def get_samsequenceline (line) -> tuple[str, int, str, str, str]:
    '''A function that splits each sequence line into seperate columns and we grab only col1 (qname), 
            col2 (bittwise flag), col3(Rname), col 4(POS), col6(CIGAR)'''
    qname = line.split("\t")[0] # gets col1 to extract UMI
    bitwiseflag = int(line.split("\t")[1]) # gets col2 for bitwiseflag
    chromname = line.split("\t")[2] # gets col3 Rname
    startpos =  line.split("\t")[3]# gets actual start position from the document
    cigar = line.split("\t")[5]# gets cigar value
    return (qname, bitwiseflag, chromname, startpos, cigar)

def get_umi(qname) -> "str":
    '''gets umi from qname which is usually at the end'''
    umi = qname.split(":")[7]
    return umi

def get_strandedness(bitwiseflag):
    '''Checks for strandedness of the read i.e is it left or right (plus or minus)'''
    if((bitwiseflag & 16) == 16): 
        strand = "-" #reverse or negative strand
    else: 
        strand = "+" #forward or positive strand
    return strand

def get_adjusted_pos(startpos, cigarstring, strand):
    '''Takes Cigar string, start position and strand and finally adjusts the start position'''
    #Reversestrand
    if strand == "-":
        onlyMDN = re.findall("(\d+)[MND]", cigarstring) #looks at digit before Match, deletion and N in cigarstring
        onlyS = re.findall("(\d+)S$", cigarstring) #looks at soft clipping at the end '$' this is important 
        sumMDN = 0
        sumS = 0
        for num in onlyMDN: #getting number from the list onlyMDN
            sumMDN += int(num)
        for number in onlyS: #getting number from the list onlyS
            sumS += int(number)
        matchingcigar = sumS + sumMDN
        adjusted_pos = int(startpos) + matchingcigar #adding into start position
    #Forwardstrand
    else:
        onlyS = re.findall("^(\d+)S", cigarstring) #we only need soft clip in the beginning denoted by ^
        sumS = 0
        for number in onlyS:
            sumS += int(number)
        adjusted_pos = int(startpos) - sumS #subtracting soft clipping from start position
    return adjusted_pos

########
# Main #
########

with open(infile,"r") as inputfile, open(umifile, 'r') as umif, open(outfile,"w") as outputfile:
    
    umi_set = set() #empty set with UMIs
    checker_set = set() #set with - umi, adjustedpos, chrm, strand
    duplicate_counter = 0 #to get the number of duplicates
    
    for umi in umif: #getting a set of umi's from STL96.txt file
        umi=umi.strip("\n")
        umi_set.add(umi) 
    #print(len(umi_set))

    for line in inputfile: #looping through each line in the file 
        if line.startswith ("@"):
            outputfile.write(line)
        else:
            qname, bitwiseflag, chromname, startpos, cigar = get_samsequenceline(line) #applying function to line
            #print(get_samsequenceline(line))
            umi = get_umi(qname) #always save the return values into a variable for function to work 
            if umi not in umi_set:
                pass 
            else:
                strand = get_strandedness(bitwiseflag) #getting strandedness
                adjusted_pos = get_adjusted_pos(startpos, cigar, strand) #getting adjusted_position
                content_checkerset = (chromname, strand, adjusted_pos, umi) #storing all the values needed in content_checkerset
                if content_checkerset in checker_set:
                    #duplicate_counter += 1
                    continue
                else:
                    checker_set.add(content_checkerset) # two parenthesis since it's a tuple
                    outputfile.write(f'{line}')
    print(duplicate_counter)        
# On command line: 
# /usr/bin/time -v ./nahata_deduper.py -f C1_SE_uniqAlign.sam -o C1_SE_uniqAlign_output.sam -u STL96.txt                    

# To look at unique reads: grep -v "^@" C1_SE_uniqAlign_output.sam | wc -l >>> 13719048
# To look at header in the output.sam file: grep "^@" C1_SE_uniqAlign_output.sam | wc -l >>> 64
# To look at unique reads per chromosome: grep -v "@" C1_SE_uniqAlign_output.sam | cut -f3 | sort -n | uniq -c > summary.txt