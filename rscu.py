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

        elif len(seq) / 100 >= 1:
            return seq
        elif len(seq) / 200 >= 1:
            return seq
        elif len(seq) / 300 >= 1:
            return seq
        elif len(seq) / 400 >= 1:
            return seq
        else:
             return None
 #codon_matrix=table.similarity_martix_codon(common_varies,codon_value_dictonary)

    def similarity_martix_codon(self,new_seq,comparer):
        #
        for i in range(0,len(new_seq),3):
            for keys,values in comparer.items():
                if new_seq[i:i+3] == keys:
                    comparer[keys] += 1
        return comparer 
 #  bablyon_tower= table.similarity_matrix_amino(babel_tower,translation_dictionary,)
    def translation(self,new_seq,translation_dict):
        proteins=''
        for i in range(0,len(new_seq),3):
                codon = new_seq[i:i+3]
                proteins += translation_dict[codon]
        return proteins

    # table.dipetide_freq(babylon_tower, dipeitide)
    def dipetide_freq(self, new_protein, dipetides):
        all_count = {}
        for dipetide in dipetides:
            count = new_protein.count(dipetide)
            if count > 0:
                all_count[dipetide] = count
        return all_count
   # table.similarity_matrix_protien(babylon_tower,protein_list):
    def similarity_matrix_protien(self,new_protein, proteins):
        all_count= {}
        for protein in proteins:
            count = new_protein.count(protein)
            if count > 0:
                all_count[protein]= count
        return all_count


    def dinucleotides_freq(self,new_seq,dinucleotides):
        all_count= {}
        for dinucleotide in dinucleotides:
            count = new_seq.count(dinucleotide)
            if count > 0:
                all_count[dinucleotide] = count
        return all_count

     def gravy_Score(self, new_protein, hydro):
        for acid in new_protein:
           GRAVYscore += hydro[acid]
         return GRAVYscore



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
    dinucleotides = ['AA', 'AT', 'AG', 'AC',
                   'TA' , 'TT', 'TG', 'TC',
                   'GA', 'GT',  'GG', 'GC',
                   'CA', 'CT', 'CG', 'CT']

    translation_dictionary= {
            'GCU': 'A' , 'GCC' : 'A' , 'GCA' : 'A' , 'GCG' : 'A' ,
            'CGU' : 'R' , 'CGC' : 'R' , 'CGA' : 'R' , 'CGG' : 'R' , 'AGA' : 'R' , 'AGG' : 'R' ,
            'UGU' : 'C' , 'UGC' : 'C',
            'AAU' : 'N' , 'AAC' : 'N' ,
            'GAU' : 'D' , 'GAC' : 'D' ,
            'GAA' : 'E' , 'GAG' : 'E' ,
            'CAA' : 'Q' , 'CAG' : 'Q' ,
            'GGU' : 'G' , 'GGC' : 'G' , 'GGA' : 'G' , 'GGG' : 'G' ,
            'CAU' : 'H' , 'CAC' : 'H' ,
            'AUU' : "I" , 'AUC' :'I' , 'AUA' : 'I' , 
            'CUU' : 'L' , 'CUC' : 'L' , 'CUA' : 'L' , 'CUG' : 'L' , 'UUA' : 'L' , 'UUG' : 'L' ,
            'AAA' : 'K' , 'AAG' : 'K' ,
            'AUG' : 'M' ,
            'UUU' : 'F' , 'UUC' : 'F' ,
            'CCU' : 'P' , 'CCC' : 'P' , 'CCA' :' P' , 'CCG' : 'P',
            'UCU' : 'S' , 'UCC' : 'S' , 'UCA' : 'S' , 'UCG' : 'S', 'AGU' : 'S' , 'AGC' : 'S' ,
            'ACU' : 'T' , 'ACC' : 'T' , 'ACA' : 'T' , 'ACG' : 'T' ,
            'UGG' : 'W' ,
            'UAU' : 'Y' , 'UAC' : 'Y' ,
            'GUU' : 'V' , 'GUC' : 'V' , 'GUA' : 'V' , 'GUG' : 'V', 
            'UAA' : ' ', 'UGA' : ' ' , 'UAG' : ' ' }
    amino_acid_list = [
            'A' , 'R' , 'N' , 'D' , 'C' ,
            'Q' , 'E' , 'G' , 'H' , 'I' ,
            'L' , 'K' , 'M' , 'F' , 'P' ,
            'S' , 'T' , 'W' , 'Y' , 'V' ,
            " " 
            ]
    dipeitide= [ 
              'AA','AR','AN','AD','AE','AQ','AG','AH','AI','AL','AK','AM','AL','AK','AF','AP','AS','AT','AW','AY','AV',"A ",
              'RA' , 'RR' , 'RC', 'RN', 'RD', 'RE' , 'RQ' , 'RG' , 'RH' , 'RI' , 'RM', 'RL' , 'RK' , 'RM' , 'RF' , 'RP' , 'RS' ,'RT', 'RW' , 'RY' , 'RV', "R ",
              'CA' , 'CR' , 'CC' , 'CN', 'CD', 'CE' , 'CQ' ,'CG', 'CH' ,'CI','CM' ,'CL', 'CK', 'CM', 'CF', 'CP' , 'CS', 'CT', 'CW', 'CY', 'CV', 'C ',
              'NA', 'NR', 'NC', 'ND' , 'NE' , 'NQ' , 'NG' ,'NH', 'NI','NM', 'NL', 'NK' , 'NN' ,'NF', 'NP','NS','NT','NW','NY','NV', 'N ',
              'DA' , 'DR' , 'DC' , 'DN' , 'DD' , 'DE' , 'DQ' , 'DG' , 'DH' , 'DI' ,  'DL' , 'DM' , 'DK', 'DF' , 'DP' , 'DS' , 'DT' , 'DW' , 'DY' , 'DV' ,'D ',
              'EA' , 'ER' , 'EC' , 'EN' , 'ED' , 'EE' , 'EQ' , 'EG' , 'EH' , 'EI' ,  'EL' , 'EM' , 'EK', 'EF' , 'EP' , 'ES' , 'ET' , 'EW' , 'EY' , 'EV' ,'E ',
              'QA' , 'QR' , 'QC' , 'QN' , 'QD' , 'QE' , 'QQ' , 'QG' , 'QH' , 'QI' ,  'QL' , 'QM' , 'QK', 'QF' , 'QP' , 'QS' , 'QT' , 'QW' , 'QY' , 'QV' ,'Q ',
              'GA' , 'GR' , 'GC' , 'GN' , 'GD' , 'GE' , 'GQ' , 'GG' , 'GH' , 'GI' ,  'GL' , 'GM' , 'GK', 'GF' , 'GP' , 'GS' , 'GT' , 'GW' , 'GY' , 'GV' ,'G ',
              'HA' , 'HR' , 'HC' , 'HN' , 'HD' , 'HE' , 'HQ' , 'HG' , 'HH' , 'HI' ,  'HL' , 'HM' ,'HK', 'HF' , 'HP' , 'HS' , 'HT' , 'HW' , 'HY' , 'HV' ,'H ',
              'IA' , 'IR' , 'IC' , 'IN' , 'ID' , 'IE' , 'IQ' , 'IG' , 'IH' , 'II' ,  'IL' , 'IM', 'IK', 'IF' , 'IP' , 'IS' , 'IT' , 'IW' , 'IY' , 'IV' ,'I ',
              'FA' , 'FR' , 'FC' , 'FN' , 'FD' , 'FE' , 'FQ' , 'FG' , 'FH' , 'FI' ,  'FL' , 'FM', 'FK', 'FF' , 'FP' , 'FS' , 'FT' , 'FW' , 'FY' , 'FV' ,'F ',
              'MA' , 'MR' , 'MC' , 'MN' , 'MD' , 'ME' , 'MQ' , 'MG' , 'MH' , 'MI' ,  'ML' , 'MM', 'MK', 'MF' , 'MP' , 'MS' , 'MT' , 'MW' , 'MY' , 'MV' ,'M ',
              'PA' , 'PR' , 'PC' , 'PN' , 'PD' , 'PE' , 'PQ' , 'PG' , 'PH' , 'PI' ,  'PL' , 'PM', 'PK', 'PF' , 'PP' , 'PS' , 'PT' , 'PW' , 'PY' , 'PV' ,'P ',
              'YA' , 'YR' , 'YC' , 'YN' , 'YD' , 'YE' , 'YQ' , 'YG' , 'YH' , 'YI' ,  'YL' , 'YM', 'YK', 'YF' , 'YP' , 'YS' , 'YT' , 'YW' , 'YY' , 'YV' ,'Y ',
              ' A' , ' R' , ' C' , ' N' , ' D' , ' E' , ' Q' , ' G' , ' H' , ' I' ,  ' L' , ' M', ' K', ' F' , ' P' , ' S' , ' T' , ' W' , ' Y' , ' V' ,'  ',
              'VA' , 'VR' , 'VC' , 'VN' , 'VD' , 'VE' , 'VQ' , 'VG' , 'VH' , 'VI' ,  'VL' , 'VM', 'VK', 'VF' , 'VP' , 'VS' , 'VT' , 'VW' , 'VY' , 'VV' ,'V ',
              ]

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
        "UAG" : 0.0, "UAA" : 0.0, "UGA" : 0.0
        }
    seq='ATGAATCCAAATCAAAAAATAATAACCATTGGATCAATCAGTATAGCAATCGGGATAATTAGTCTAATGTTGCAAATAGGAAATATTATTTCAATATGGGCTAGTCACTCAATCCAAACTGGAAGTCAAAACCACACTGGAATATGCAACCAAAGAATCATCACATATGAAAACAGCACCTGGGTGAATCACACATATGTTAGTATTAACAACACTAATGTTGTTGCTGGAAAGGACAAAACTTCAGTGACATTGGCCGGCAATTCATCTCTTTGTTCTATCAGTGGATGGGCTATATACACAAAAGACAACAGCATAAGAATTGGCTCCAAAGGAGATGTTTTTGTCATAAGAGAACCTTTCATATCATGTTCTCATTTGGAATGCAGGACCTTTTTTCTGACCCAAGGTGCTCTATTAAATGACAAACATTCAAATGGAACCGTTAAGGACCGAAGTCCTTATAGGGCCCTAATGAGCTGTCCTCTAGGTGAAGCTCCGTCTCCATACAATTCAAAGTTTGAATCAGTTGCATGGTCAGCAAGCGCATGCCATGATGGCATGGGCTGGTTAACAATCGGAATTTCTGGTCCAGACAATGGAGCTGTGGCTGTACTAAAATACAACGGAATAATAACTGAAACCATAAAAAGTTGGAAAAAGCGAATATTGAGAACACAAGAGTCTGAATGTGTCTGTGTGAACGGGTCATGTTTCACCATAATGACCGATGGCCCGAGTAATGGGGCCGCCTCGTACAAAATCTTCAAGATCGAAAAGGGGAAGGTTACTAAATCAATAGAGTTAAATGCACCCAATTTTCATTATGAGGAATGTTCCTGTTACCCAGACACTGGCACAGTGATGTGTGTATGCAGGGACAACTGGCATGGTTCAAATCGACCTTGGGTGTCTTTTAATCAAAACTTGGATTATCAAATAGGATACATCTGCAGTGGAGTGTTCGGTGACAATCCGCGTCCCAAAGATGGAAAGGGCAGCTGTAATCCAGTGACTGTTGATGGAGCAGACGGGGTTAAGGGGTTTTCATACAAATATGGTAATGGTGTTTGGATAGGAAGAACTAAAAGTAACAGACTTAGAAAGGGGTTTGAGATGATTTGGGATCCTAATGGATGGACAGATACCGACAGTGATTTCTCAGTGAAACAGGATGTTGTGGCAATAACTGATTGGTCAGGGTACAGTGGAAGTTTCGTTCAACATCCTGAGTTAACAGGATTGGACTGTATAAGACCTTGCTTCTGGGTTGAGTTAGTCAGAGGACTGCCTAGAGAAAATACAACAATCTGGACTAGTGGGAGCAGCATTTCTTTTTGTGGCGTTGATAGTGATACTGCAAACTGGTCTTGGCCAGACGGTGCTGAGTTGCCGTTCACCATTGACAAGTAG'
    seq1='ATGAAAGTAAAACTACTGGTCCTATTATGCACATTTACAGCTACATATGCAGACACAATATGTATCGGCTACCATGCCAACAACTCAACCGACACTGTTGACACAGTACTTGAAAAGAATGTGACAGTGACACACTCTGTCAACCTGCTTGAGGACAACCACAATGGAAAACTATGTCTATTAAAAGGAATAGCCCCATTACAATTGGGTAACTGCAGCGTTGCCGGGTGGATCTTAGGAAACCCAGAATGCGAATTACTGATTTCCAAGGAGTCATGGTCCTACATTGTAGAAAAACCAAATCCTGAGAATGGAACATGTTACCCAGGGCATTTCGCCGACTATGAGGAACTGAGGGAGCAATTGAGTTCAGTATCTTCATTTGAGAGGTTCGAAATATTCCCCAAAGAAAGCTCATGGCCCAACCACACCGTAACCGGAGTATCCGCATCATGCTCCCATAATGGGGAAAGCAGCTTTTACAAAAATTTGCTATGGCTGACGGGAAAGAATGGTTTGTACCCAAACCTGAGCAAGTCCTATGCAAACAACAAAGAGAAAGAAGTCCTCGTACTATGGGGTGTTCATCACCCGCCAAACATAGGTGACCAAATGACCCTCTATCATAAAGAAAATGCTTATGTCTCTGTAGTGTCTTCACATTATAGCAGAAAATTCACCCCAGAAATAGCCAAAAGACCCAAAGTAAGAGATCAAGAAGGAAGAATCAACTACTACTGGACTCTGCTTGAACCCGGGGATACAATAATATTTGAGGCAAATGGAAATCTAATAGCGCCAAGATATGCTTTCGCACTGAGTAGAGGCTTTGGATCAGGAATCATCAACTCAAATGCACCAATGGATGAATGTGATGCGAAGTGCCAAACACCTCAGGGAGCTATAAACAGCAGTCTTCCTTTCCAGAATGTACACCCAGTCACAATAGGAGAATGTCCAAAGTATGTCAGGAGTGCAAAATTAAGGATGGTTACAGGACTAAGGAACATCCCATCCATTCAATCCAGAGGTTTGTTTGGAGCCATTGCCGGTTTCATTGAAGGGGGGTGGACTGGAATGGTAGATGGTTGGTATGGTTATCATCACCAGAATGAGCAAGGATCTGGCTATGCTGCAGATCAAAAAAGCACACAAAATGCCATTAATGGGATTACAAACAAGGTGAATTCTGTAATTGAGAAGATGAACACTCAATTCACAGCTGTGGGCAAAGAATTCAACAAATTGGAAAGAAGGATGGAAAACTTAAATAAAAAAGTTGATGACGGGTTTATAGACGTTTGGACATATAATGCAGAACTGTTGGTTCTACTGGAAAATGAAAGAACTTTGGATTTCCATGACTCCAATGTGAAGAATTTGTATGAGAAAGTAAAAAATCAATTAAAGAATAATGCCAAAGAAATAGGAAATGGGTGTTTTGAATTTTATCACAAGTGTAACGATGAATGCATGGAGAGTGTAAAAAATGGAACTTATGACTATCCAAAATATTCCGAAGAATCAAAGTTAAGCAGGGAGAAAATTGATGGAGTGAAATTGGAATCAATGGGAGTCTATCAGATTCTGGCGATCTACTCAACAGTCGCCAGTTCTCTGGTTCTTTTGGTCTCCCTGGGGGCAATCAGCTTCTGGATGTGTTCCAATGGGTCTTTGCAGTGTAGAATATGCATCTAA'
    
    hydropathy = {"A":1.800,"R":-4.500,"N":-3.500,"D":-3.500,"C":2.500,"Q":-3.500,"E":-3.500,"G":-0.400,"H":-3    .200,"I":4.500,"L":3.800,"K":-3.900,"M":1.900,"F":2.800,"P":-1.600,"S":-0.800,"T":-0.700,"W":-0.900,"Y":-1.300,"V" : 4.200} 
    table=mapping(seq1)

    common_varies=table.removing_common_comp()

    dinucleotide_freqency= table.dinucleotides_freq(common_varies,dinucleotides)

    common_varies = common_varies.replace('T','U')
    
    tabling=mapping(seq)
    not_optimal= tabling.removing_common_comp()
    not_optimal= not_optimal.replace('T','U')
    
    print ("%s \n",common_varies)
    print ("%s \n",not_optimal)

    babel=(table.divides_by_three(common_varies))
    non_optimal=(tabling.divides_by_three(not_optimal))
    print(dinucleotide_freqency) 


    print("%s \n",babel)
    print("%s \n",non_optimal)
    bablyon_tower= table.translation(common_varies,translation_dictionary)
    protein_frequncy = table.similarity_matrix_protien(bablyon_tower,amino_acid_list)

    def_optimal=table.greater_then_hundred(non_optimal)
    dipetide_frequency = table.dipetide_freq(bablyon_tower,dipeitide)
    tower=(table.greater_then_hundred(common_varies))

    print("%s \n",tower)
    print("%s \n",def_optimal)

    babel_tower=table.similarity_martix_codon(common_varies,codon_value_dictonary)

    def_optimized = tabling.similarity_martix_codon(def_optimal,codon_value_dictonary)
    
    sort_bable_tower= sorted(babel_tower.items(),key =lambda x : x[1], reverse = True)

    sort_def_optimized= sorted(def_optimized.items(),key =lambda x : x[1], reverse = True)

    for i in sort_bable_tower:
       print("%s : %.1f" %(i[0], i[1]))
    for i in sort_def_optimized:
       print("%s : %.1f"%(i[0],i[1]))

    #print(bablyon_tower)

    #print(type(protein_frequncy))

    #print(protein_frequncy)
    #sorted_dipetide_freq = sorted(dipetide_frequency.items(),key =lambda x : x[1], reverse = True)
    #for i in sorted_dipetide_freq:
    #   print ( " %s : %.1f"%(i[0],i[1]))

if __name__=='__main__':
    main()
