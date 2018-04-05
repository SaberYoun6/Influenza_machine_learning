#!/etc/bin/env python 
# Samuel Young
## RSCU.py
## 
## 
#from biopy import seq 
import pdb
##gobal varialbes #####

foo = 0x0076
bar = 0x0023
shell= 0x0034




class mapping(object):
    #def __init__(self,sequences, animo_acids,codon):
    #    self.sequences= sequences
    #    self.animo_acids= amino_acid
    #    self.codons=condons
    def  removing_common_comp(self,sequences,stopCodon=('TAA','TAG','TGA')):
        seq=''
        seqs=''
        sequences.rstrip()
        if (sequences.startswith('ATG')):
            seq=sequences.replace('ATG','')
        if (seq.endswith('TAA')):
              seqs = seq.replace('TAA','')
        elif seq.endwith('TAG'):
              seqs = seq.replace('TAG','')
        else:
              seqs = seq.replace('TGA','')
        return seqs
    def divides_by_three(self,seq):
        if(len(seq)%3 == 0):
	    return True
        else:
            return False
    #def similarity_martix(self,sequences,comparer,similar_results):
    #    pdb.set_trace()
    #    count = 0
    #    for seq in sequences:
    #        for keys,values in comparer.items():
    #            for key,value in similar_results.items():
    #                if (seq is values):
    #                if (seq is key):
    #                while True:
    #                    if (seq.index(keys)  seq.index(key)):
    #                        seq[sim]=+seq[comp]
    #                        number=len(similar_results[key])

    #def rel_syn_cod_use(self,sequences,amino_acids,codons):
    #    pdb.set_trace()
    #    j=0
    #    for a,i in animo_acids:
    #        for c,j in codons:
    #            for X in enumerate(sequences):
    #                lenOriginal=X[i][j] 
    #                print (lenOriginal)
    #                sumOriginal= sum(X[len(len(sequences)*i)][j-1])
    #                rel_syn_cod_usage=lenOriginal/((1/number)*sumOriginal)
    #    return rel_syn_cod_usage

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

def  main():
    table=mapping()
    translation_dictionary= {
            'Ala' : 'GCT' ,  'Ala': 'GCC' , 'ALA' : 'GCA' , 'Ala' : 'GCG' ,
            'Arg' : 'CGT' ,' Arg' : 'CGC' , 'Arg' : 'CGA' , 'Arg' : 'CGG' , 'Arg' : 'AGA' , 'Arg' : 'AGG' ,
            'Asn' : 'AAT' ,' Asn' : 'AAC' , 'Asp' : 'GAT' , 'Asp' : 'GAC' ,
            'Cys' : 'TGT' , 'Cys' : 'TGC' , 'Gln' : 'CAA' , 'Gln' : 'CAG' , 'Glu' : 'GAA' , 'Glu ': 'GAG' ,
            'Gly' : 'GGT' , 'Gly' : 'GGC' , 'Gly' : 'GCA' , 'Gly' : 'GCG' , 'His' : 'CAT' , 'His' : 'CAC' ,
            'Ile' : "ATT" , 'Ile' :' ATC' , 'Ile' : 'ATA' , 
            'Leu' : 'CTT' , 'Leu' : 'CTC' , 'Leu' : 'CTA' , 'Leu' : 'CTG' , 'Leu' : 'TTA' , 'Leu' : 'TTG' ,
            'Lys' : 'AAA' , 'Lys' : 'AAG' , 'Met' : 'ATG' , 'Phe' : 'TTT' ,' Phe' : 'TTC' ,
            'Pro' : 'CCT' , 'Pro' : 'CCC' , 'Pro' :' CCA' , 'Pro' : 'CCG',
            'Ser' : 'TCT' , 'Ser' : 'TCC' , 'Ser' : 'TCA' , 'Ser' : 'TCG', 'Ser' : 'AGT' , 'Ser' : 'AGC' ,
        'Thr' : 'CGT' , 'Thr' : 'ACC' , 'Thr' : 'ACA' , 'Thr' : 'ACG' , 'Trp' : 'TGG',
        'Tyr' : 'TAT' , 'Tyr' : 'TAC' , 'Val' : 'GTT' , 'Val' : 'GTC' , 'Val' : 'GTA' , 'Val' : 'GTG'
        }
    amino_acid_dictionary = {
            'Ala' : 0.0 , 'Arg' : 0.0  , 'Asn' : 0.0 , 'Asp' : 0.0 , 'Cys' : 0.0 ,
            'Gln' : 0.0 , 'Glu' : 0.0  , 'Gly' : 0.0 , 'His' : 0.0 , 'Ile' : 0.0 ,
            'Leu' : 0.0 , 'Lys' : 0.0  , 'Met' : 0.0 , 'Phe' : 0.0 , 'pro' : 0.0 ,
            'Ser' : 0.0 , 'Thr' : 0.0  , 'Trp' : 0.0 , 'Tyr' : 0.0 , 'Val' : 0.0
            }
    codon_value_dictonary ={
        'GCT' : 0.0 , 'GCC' : 0.0 , 'GCA' : 0.0 , 'GCG' : 0.0 ,
        'CGT' : 0.0 , 'CGC' : 0.0 , 'CGA' : 0.0 , 'CGG' : 0.0 , 'AGA' : 0.0 , 'AGG' : 0.0 ,
        'AAT' : 0.0 , 'AAC' : 0.0 , 'GAT' : 0.0 , 'GAC' : 0.0 ,
        'TGT' : 0.0 , 'TGC' : 0.0 , 'CAA' : 0.0 , 'CAG' : 0.0 , 'GAA' : 0.0 , 'GAG' : 0.0 ,
        'ATT' : 0.0 , 'ATC' : 0.0 , 'ATA' : 0.0 ,
        'CTT' : 0.0 , 'CTC' : 0.0 , 'CTA' : 0.0 , 'CTG' : 0.0 , 'TTA' : 0.0 , 'TTG' : 0.0 ,
        'AAA' : 0.0 , 'AAG' : 0.0 , 'ATG' : 0.0 , 'TTT' : 0.0 , 'TTC' : 0.0 ,
        'CCT' : 0.0 , 'CCC' : 0.0 , 'CCA' : 0.0 , 'CCG' : 0.0 ,
        'TCT' : 0.0 , 'TCC' : 0.0 , 'TCA' : 0.0 , 'TCG' : 0.0 , 'AGT' : 0.0 , 'AGC' : 0.0 ,
        'CGT' : 0.0 , 'ACC' : 0.0 , 'ACA' : 0.0 , 'ACG' : 0.0,  'TGG' : 0.0 ,
        'TAT' : 0.0 , 'TAC' : 0.0 , 'GTT' : 0.0 , 'GTC' : 0.0,  'GTA' : 0.0 , 'GTG' : 0.0 ,
        }
    seq='ATGCTGTCTCAGGCACGTGGATGGTTTGGACAAATCAGATTCAAGTCTGATCAACCTTCATACAGATCTAGAGTCTAAAGCAGTTAA' 
    seq1='AATTTTCCCACTGCCTTAAGCCGGCTTGCCCTTTCTGCCTGTAGATCCATTGGACTGGTGCCAACGCGCAGGCATAGTTCGAGGAGAATTATCCGGGGGCAATGACAACCAGCATCTCGGGTCTTGCCCAACCCGTCTACACGCTGTTATAGCGTATCAGCGGGAACCCGGTGCCACGCGATGGAACGTCCAACTCTGGCAGGCAATTAAAGGGAACGTA'
    
    common_varies=table.removing_common_comp(seq)
    seq1=seq1[::-1]
    not_optimal= table.removing_common_comp(seq1)
    print (common_varies)
    print (not_optimal)
    babel=table.divides_by_three(common_varies)
    non_optimal=table.divides_by_three(not_optimal)
    #babel=table.similarity_martix(seq,translation_dictionary,amino_acid_dictionary)
    print babel
    print non_optimal 


if __name__=='__main__':
    main()
