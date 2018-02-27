#import dna_calculatortests as dna

import itertools
import protein_calculatorstest as prt
import sys

class method():
    #def __init__(self):

    def read_file(self,filename):
    # this should allow us to read in a file and output the result that we want
        head_seqs={}
        with open(filename) as f:
            #we are reading in the file
            headerfile=f.next().strip()
            headerfile=headerfile.split('\n')
            seq=f.next().strip()
            seq=seq.split('\n')
            # this should return all the value in head_seq{}
            for k,v in itertools.groupby(headerfile,lambda x: x[len(seq)]):
                    head_seqs.setdefault(k,[]).append(v)
            return head_seqs








def main():
    mthd= method()
    name =sys.argv[1]
    seqdata=mthd.read_file(name)
    print seqdata
 #   dna.fullProfile(seqdata)

main()
