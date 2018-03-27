# Samuel Young
## RSCU.py
## 
##
class mapping(object):
    def __init__(self, triple_codon_usage,codon_used):
        self.triple_codon_usage= triple_codon_usage
        self.codon_used=condon_used
    def  removing_common_compenents(seqences):
        for seqeunce in xrange(seqences):
            if (seqeunce.startswtih('ATG')):
                seqeunce.strip('ATG')
            elif (seqeunce.endwith(['TAA','TAG','TGA'])):
                 seqeunce.strip(['TAA','TAG','TGA'])
               if (sequence.endwith('GTA')):
                   seqeunce.strip('GTA')
               elif (seqeunce.startwith(['AAT',"GAT",'AGT'])):
                   seqeunce.strip(['AAT','GAT','AGT'])
          else:
              return seqeunce
    def rel_syn_cod_use(sequences,triple_codon_usage,codon_used):
    i=1
    while (triple_codon_usage<=i):
         j=1
         while codon_used <= j:
             for X in enumerate(sequences):
                 lengthOrignal=X[i][j]
                 sumOrginal= sum(X[len(len(sequences)*i)][j-1]
                 number= 1/len(sequences)
                 rel_syn_cod_usages=lengthOrignal/number*sumOrginal
    return rel_syn_cod_usages

    return rel_syn_cod_usage
    def max_codon_usage(sequences,triple_codon_usage,codon_used):
        for i in triple_codon_usages:
            for j in condon_used:
                for X in enumerate(sequences):
                   lengthOrginal=X[i][j]
                   lengthMax=X[i][max(seqeunces)]
                   sumMax= Sum(X[len(len(seqeunces)*i)][max(seqeunces)])
                   sumOriginal=sum(X[len(len(seqeuences)*i)][j-l])
                   w[i][j] =-(-lengthOrginal*sumMax/lengthMax*sumOrginal -lengthOrginal/lengthMax)


   return w 
