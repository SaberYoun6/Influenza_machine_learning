#!/etc/bin/env python 
# Samuel Young
## RSCU.py
## 
## 
#from biopy import seq 

##gobal varialbes #####

foo = 0x76
bar = 0x23
shell= 0x34




class mapping(object):
    #def __init__(self, triple_codon_usage,codon_used):
    #    self.triple_codon_usage= triple_codon_usage
    #    self.codon_used=condon_used
    def  removing_common_comp(self,seqences,stopCodon=('TAA','TAG','TGA')):
        
        for seqeunce in seqences:
            if (seqeunce.startswith('ATG')):
                seqeunce.strip('ATG')
                print ('what going on?')
            elif (seqeunce.endswith(stopCodon)):
                if (seqeunce.endwith('TAA')):
                    seqeunce.index('TAA')
                    seqeunce.strip('TAA')
                    print (foo)
                elif (sequnce.endwith('TAG')):
                    seqeunce.index('TAG')
                    seqeune.strip('TAG')
                    print(bar)
                else:
                    seqeunce.index('TGA')
                    seqeunce.strip('TGA')
                    print(shell)
            

    
        return seqeunce
    #''''
    #def rel_syn_cod_use(sequences,triple_codon_usage,codon_used):
    #    i=0
    #    while (len(sequences)/3<=i):
    #        j=0
    #        while (len(sequences) <= j):
    #            number= len(sequences)
    #            for X in enumerate(sequences):
    #                lengthOriginal=X[i][j]
    #                j=+1
    #                i=+1
    #                print (lengthOriginal)
    #                sumOriginal= sum(X[len(len(sequences)*i)][j-1]
    #                        rel_syn_cod_usages=lengthOriginal/((1/number)*sumOriginal)
    #    return rel_syn_cod_usages

    #def max_codon_usage(sequences,triple_codon_usage,codon_used):
    #    for i in triple_codon_usages:
    #        for j in condon_used:
    #            for X in enumerate(sequences):
    #               lengthOrginal=X[i][j]
    #               lengthMax=X[i][(max(len(seqeunces)/3)]
    #               sumMax= Sum(X[len(len(seqeunces)*i)][max(len(seqeunces)/3)])
    #               sumOriginal=sum(X[len(len(seqeuences)*i)][j-l])
    #               w[i][j] =-(-lengthOrginal*sumMax/lengthMax*sumOrginal -lengthOrginal/lengthMax)
    #    return w 
    # ''''

def  main():
    table=mapping()
    seq='ATGCTGTCTCAGGCACGTGGATGGTTTGGACAAATCAGATTCAAGTCTGATCAACCTTCATACAGATCTAGAGTCTAAAGCAGTTAA' 
    seq1='AATTTTCCCACTGCCTTAAGCCGGCTTGCCCTTTCTGCCTGTAGATCCATTGGACTGGTGCCAACGCGCAGGCATAGTTCGAGGAGAATTATCCGGGGGCAATGACAACCAGCATCTCGGGTCTTGCCCAACCCGTCTACACGCTGTTATAGCGTATCAGCGGGAACCCGGTGCCACGCGATGGAACGTCCAACTCTGGCAGGCAATTAAAGGGAACGTA'
    common_varies =table.removing_common_comp(seq)
    seq1=seq1[::-1]
    not_optimal= table.removing_common_comp(seq1)
    print (common_varies)
    print (not_optimal)


if __name__=='__main__':
    main()
