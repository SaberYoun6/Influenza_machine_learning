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
    def __init__(self,sequences):
       self.sequences= sequences
    def  removing_common_comp(self,stopCodon=('TAA','TAG','TGA')):
        sequences =self.sequences
        seq=''
        seqs=''
        sequences.rstrip()
        if sequences.startswith('ATG'):
            seq=sequences.replace('ATG','')
        if seq.endswith('TAA'):
              seqs = seq.replace('TAA','')
        elif seq.endswith('TAG'):
              seqs = seq.replace('TAG','')
        else:
              seqs = seq.replace('TGA','')
        return seqs

    def divides_by_three(self,new_seq):
        seq = new_seq
        if(len(seq)%3 == 0):
            return seq
        else:
            return None

    def greater_then_hundred(self,new_seq):
        seq = new_seq
        if (seq == None):
            return None

        elif(len(seq)  / 100 >= 1) :
            return seq

        else:
             return None

    def similarity_martix_codon(self,new_seq,comparer):
        
        for i in range(0,len(new_seq),3):
            for keys,values in comparer.items():
                if new_seq[i:i+3] == keys:
                    comparer[keys] += 1
        return comparer 

    def similarity_matrix_amino(self,codon,trans):
        temp_dict={}
        
        for k,v in codon.items():
            values= 0
            for ke,va in trans.items():
                if va == k :
                    values += v
                temp_dict[ke]= values

        return temp_dict 


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
    translation_dictionary= {
            'A' : 'GCU' , 'A' : 'GCC' , 'A' : 'GCA' , 'A' : 'GCG' ,
            'R' : 'CGU' , 'R' : 'CGC' , 'R' : 'CGA' , 'R' : 'CGG' , 'R' : 'AGA' , 'R' : 'AGG' ,
            'C' : 'UGU' , 'C' : 'UGC',
            'N' : 'AAU' , 'N' : 'AAC' , 'D' : 'GAU' , 'D' : 'GAC' ,
            'E' : 'GAA' , 'E' : 'GAG' , 'Q' : 'CAA' , 'Q' : 'CAG' ,
            'G' : 'GGU' , 'G' : 'GGC' , 'G' : 'GGA' , 'G' : 'GGG' , 'H' : 'CAU' , 'H' : 'CAC' ,
            'I' : "AUU" , 'I' :' AUC' , 'I' : 'AUA' , 
            'L' : 'CUU' , 'L' : 'CUC' , 'L' : 'CUA' , 'L' : 'CUG' , 'L' : 'UUA' , 'L' : 'UUG' ,
            'K' : 'AAA' , 'K' : 'AAG' , 'M' : 'AUG' , 'F' : 'UUU' , 'F' : 'UUC' ,
            'P' : 'CCU' , 'P' : 'CCC' , 'P' :' CCA' , 'P' : 'CCG',
            'S' : 'UCU' , 'S' : 'UCC' , 'S' : 'UCA' , 'S' : 'UCG', 'S' : 'AGU' , 'S' : 'AGC' ,
            'T' : 'ACU' , 'T' : 'ACC' , 'T' : 'ACA' , 'T' : 'ACG' , 'W' : 'UGG',
            'Y' : 'UAU' , 'Y' : 'UAC' , 'V' : 'GUU' , 'V' : 'GUC' , 'V' : 'GUA' , 'V' : 'GUG'
        }
    amino_acid_dictionary = {
            'A' : 0.0 , 'R' : 0.0  , 'N' : 0.0 , 'D' : 0.0 , 'C' : 0.0 ,
            'Q' : 0.0 , 'E' : 0.0  , 'G' : 0.0 , 'H' : 0.0 , 'I' : 0.0 ,
            'L' : 0.0 , 'K' : 0.0  , 'M' : 0.0 , 'F' : 0.0 , 'P' : 0.0 ,
            'S' : 0.0 , 'T' : 0.0  , 'W' : 0.0 , 'Y' : 0.0 , 'V' : 0.0
            }
    codon_value_dictonary ={
        'GCU' : 0.0 , 'GCC' : 0.0 , 'GCA' : 0.0 , 'GCG' : 0.0 ,
        'CGU' : 0.0 , 'CGC' : 0.0 , 'CGA' : 0.0 , 'CGG' : 0.0 , 'AGA' : 0.0 , 'AGG' : 0.0 ,
        'AAU' : 0.0 , 'AAC' : 0.0 , 'GAU' : 0.0 , 'GAC' : 0.0 ,
        'UGU' : 0.0 , 'UGC' : 0.0 , 'CAA' : 0.0 , 'CAG' : 0.0 , 'GAA' : 0.0 , 'GAG' : 0.0 ,
        'AUU' : 0.0 , 'AUC' : 0.0 , 'AUA' : 0.0 ,
        'CUU' : 0.0 , 'CUC' : 0.0 , 'CUA' : 0.0 , 'CUG' : 0.0 , 'UUA' : 0.0 , 'UUG' : 0.0 ,
        'AAA' : 0.0 , 'AAG' : 0.0 , 'AUG' : 0.0 , 'UUU' : 0.0 , 'UUC' : 0.0 ,
        'CCU' : 0.0 , 'CCC' : 0.0 , 'CCA' : 0.0 , 'CCG' : 0.0 ,
        'UCU' : 0.0 , 'UCC' : 0.0 , 'UCA' : 0.0 , 'UCG' : 0.0 , 'AGU' : 0.0 , 'AGC' : 0.0 ,
        'CGU' : 0.0 , 'ACC' : 0.0 , 'ACA' : 0.0 , 'ACG' : 0.0,  'UGG' : 0.0 ,
        'UAU' : 0.0 , 'UAC' : 0.0 , 'GUU' : 0.0 , 'GUC' : 0.0,  'GUA' : 0.0 , 'GUG' : 0.0 ,
        }
    seq='ATGAATCCAAATCAAAAAATAATAACCATTGGATCAATCAGTATAGCAATCGGGATAATTAGTCTAATGTTGCAAATAGGAAATATTATTTCAATATGGGCTAGTCACTCAATCCAAACTGGAAGTCAAAACCACACTGGAATATGCAACCAAAGAATCATCACATATGAAAACAGCACCTGGGTGAATCACACATATGTTAGTATTAACAACACTAATGTTGTTGCTGGAAAGGACAAAACTTCAGTGACATTGGCCGGCAATTCATCTCTTTGTTCTATCAGTGGATGGGCTATATACACAAAAGACAACAGCATAAGAATTGGCTCCAAAGGAGATGTTTTTGTCATAAGAGAACCTTTCATATCATGTTCTCATTTGGAATGCAGGACCTTTTTTCTGACCCAAGGTGCTCTATTAAATGACAAACATTCAAATGGAACCGTTAAGGACCGAAGTCCTTATAGGGCCCTAATGAGCTGTCCTCTAGGTGAAGCTCCGTCTCCATACAATTCAAAGTTTGAATCAGTTGCATGGTCAGCAAGCGCATGCCATGATGGCATGGGCTGGTTAACAATCGGAATTTCTGGTCCAGACAATGGAGCTGTGGCTGTACTAAAATACAACGGAATAATAACTGAAACCATAAAAAGTTGGAAAAAGCGAATATTGAGAACACAAGAGTCTGAATGTGTCTGTGTGAACGGGTCATGTTTCACCATAATGACCGATGGCCCGAGTAATGGGGCCGCCTCGTACAAAATCTTCAAGATCGAAAAGGGGAAGGTTACTAAATCAATAGAGTTAAATGCACCCAATTTTCATTATGAGGAATGTTCCTGTTACCCAGACACTGGCACAGTGATGTGTGTATGCAGGGACAACTGGCATGGTTCAAATCGACCTTGGGTGTCTTTTAATCAAAACTTGGATTATCAAATAGGATACATCTGCAGTGGAGTGTTCGGTGACAATCCGCGTCCCAAAGATGGAAAGGGCAGCTGTAATCCAGTGACTGTTGATGGAGCAGACGGGGTTAAGGGGTTTTCATACAAATATGGTAATGGTGTTTGGATAGGAAGAACTAAAAGTAACAGACTTAGAAAGGGGTTTGAGATGATTTGGGATCCTAATGGATGGACAGATACCGACAGTGATTTCTCAGTGAAACAGGATGTTGTGGCAATAACTGATTGGTCAGGGTACAGTGGAAGTTTCGTTCAACATCCTGAGTTAACAGGATTGGACTGTATAAGACCTTGCTTCTGGGTTGAGTTAGTCAGAGGACTGCCTAGAGAAAATACAACAATCTGGACTAGTGGGAGCAGCATTTCTTTTTGTGGCGTTGATAGTGATACTGCAAACTGGTCTTGGCCAGACGGTGCTGAGTTGCCGTTCACCATTGACAAGTAG'
    seq1='ATGAAAGTAAAACTACTGGTCCTATTATGCACATTTACAGCTACATATGCAGACACAATATGTATCGGCTACCATGCCAACAACTCAACCGACACTGTTGACACAGTACTTGAAAAGAATGTGACAGTGACACACTCTGTCAACCTGCTTGAGGACAACCACAATGGAAAACTATGTCTATTAAAAGGAATAGCCCCATTACAATTGGGTAACTGCAGCGTTGCCGGGTGGATCTTAGGAAACCCAGAATGCGAATTACTGATTTCCAAGGAGTCATGGTCCTACATTGTAGAAAAACCAAATCCTGAGAATGGAACATGTTACCCAGGGCATTTCGCCGACTATGAGGAACTGAGGGAGCAATTGAGTTCAGTATCTTCATTTGAGAGGTTCGAAATATTCCCCAAAGAAAGCTCATGGCCCAACCACACCGTAACCGGAGTATCCGCATCATGCTCCCATAATGGGGAAAGCAGCTTTTACAAAAATTTGCTATGGCTGACGGGAAAGAATGGTTTGTACCCAAACCTGAGCAAGTCCTATGCAAACAACAAAGAGAAAGAAGTCCTCGTACTATGGGGTGTTCATCACCCGCCAAACATAGGTGACCAAATGACCCTCTATCATAAAGAAAATGCTTATGTCTCTGTAGTGTCTTCACATTATAGCAGAAAATTCACCCCAGAAATAGCCAAAAGACCCAAAGTAAGAGATCAAGAAGGAAGAATCAACTACTACTGGACTCTGCTTGAACCCGGGGATACAATAATATTTGAGGCAAATGGAAATCTAATAGCGCCAAGATATGCTTTCGCACTGAGTAGAGGCTTTGGATCAGGAATCATCAACTCAAATGCACCAATGGATGAATGTGATGCGAAGTGCCAAACACCTCAGGGAGCTATAAACAGCAGTCTTCCTTTCCAGAATGTACACCCAGTCACAATAGGAGAATGTCCAAAGTATGTCAGGAGTGCAAAATTAAGGATGGTTACAGGACTAAGGAACATCCCATCCATTCAATCCAGAGGTTTGTTTGGAGCCATTGCCGGTTTCATTGAAGGGGGGTGGACTGGAATGGTAGATGGTTGGTATGGTTATCATCACCAGAATGAGCAAGGATCTGGCTATGCTGCAGATCAAAAAAGCACACAAAATGCCATTAATGGGATTACAAACAAGGTGAATTCTGTAATTGAGAAGATGAACACTCAATTCACAGCTGTGGGCAAAGAATTCAACAAATTGGAAAGAAGGATGGAAAACTTAAATAAAAAAGTTGATGACGGGTTTATAGACGTTTGGACATATAATGCAGAACTGTTGGTTCTACTGGAAAATGAAAGAACTTTGGATTTCCATGACTCCAATGTGAAGAATTTGTATGAGAAAGTAAAAAATCAATTAAAGAATAATGCCAAAGAAATAGGAAATGGGTGTTTTGAATTTTATCACAAGTGTAACGATGAATGCATGGAGAGTGTAAAAAATGGAACTTATGACTATCCAAAATATTCCGAAGAATCAAAGTTAAGCAGGGAGAAAATTGATGGAGTGAAATTGGAATCAATGGGAGTCTATCAGATTCTGGCGATCTACTCAACAGTCGCCAGTTCTCTGGTTCTTTTGGTCTCCCTGGGGGCAATCAGCTTCTGGATGTGTTCCAATGGGTCTTTGCAGTGTAGAATATGCATCTAA'
    
    table=mapping(seq)
    common_varies=table.removing_common_comp()
    common_varies = common_varies.replace('T','U')
    
    tabling=mapping(seq1)
    not_optimal= tabling.removing_common_comp()
    not_optimal= not_optimal.replace('T','U')
    
    print ("%s \n",common_varies)
    print ("%s \n",not_optimal)


    babel=(table.divides_by_three(common_varies))
    non_optimal=(tabling.divides_by_three(not_optimal))
    
    babel_tower=table.similarity_martix_codon(common_varies,codon_value_dictonary)




    print("%s \n",babel)
    print("%s \n",non_optimal)
    bablyon_tower= table.similarity_matrix_amino(babel_tower,translation_dictionary)

    tower=(table.greater_then_hundred(common_varies))
    def_optimal=table.greater_then_hundred(non_optimal)

    print("%s \n",tower)
    print("%s \n",def_optimal)

    print(babel_tower)
    
    print(bablyon_tower)

if __name__=='__main__':
    main()
