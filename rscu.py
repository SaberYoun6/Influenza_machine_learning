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
            rel_syn_cod_usage[i][j]= X[i][j]/(1/len(sequences)*sum(X[(len(seq)*i))][(j-1)])

    return rel_syn_cod_usage
    def max_codon_usage(sequences,triple_codon_usage,codon_used):
        for i in triple_codon_usages:
            for j in condon_used:
                for X in enumerate(sequences):
                    w[i][j]= X[i][j] = x[i][j]* sum(X[len(sequences)*i)][max(seqeunces/3)-1)]/(X[i][max(sequences/3)*(sum(X[

