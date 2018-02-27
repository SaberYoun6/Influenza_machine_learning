#import dna_calculatortests as dna

import itertools
#import protein_calculatorstest as prt
import sys
### variable ####
#No Global
###


class method():
    #def __init__(self):

    def read_file(self,filename):
    # this should allow us to read in a file and output the result that we want
        ####variables
        head_seqs={}
        
        ###
        with open(filename) as f:
            #we are reading in the file
#the next few lines are used to obtain both the header of the file and then the seq of the file as well 
            fh= f.readline()
            headerfile=fh.split("\n")
            if(headerfile().startswith('>')):
                header = headerfile 
            else:
               headerfile =seq
            for header, seq in enumerate(f):
                head_seqs['header'+str(header)]= [seq] 
            return head_seqs








def main():
    mthd= method()
    name =sys.argv[1]
    seqdata=mthd.read_file(name)
    print seqdata
    #        print(key + ","+ value + "\n")
    #dna.fullProfile(seqdata)

main()
