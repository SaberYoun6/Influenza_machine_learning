import DNA_calculators as DNA
import numpy
import math
import sys

inputSequences = sys.argv[1] 
outputAnalysis = sys.argv[2] 

fsequences = open(inputSequences,"r")
sequencesData = fsequences.readlines()
fsequences.close()


headers = []
sequences = []
gene = ""
i=0
for line in sequencesData:
	if ">" in line:
		headers.append(line.strip())
		if i>0:
			sequences.append(gene)
			gene = ""
	else:
		gene += line.strip()
	i+=1
sequences.append(gene)
	

HAheaders = []
HAsequencesDNA = []
NAheaders = []
NAsequencesDNA = []

for header in headers:
	if "Segment:4" in header:
		HAheaders.append(header)
		index = headers.index(header)
		HAsequencesDNA.append(sequences[index])

	if "Segment:6" in header:
		NAheaders.append(header)
		index = headers.index(header)
		NAsequencesDNA.append(sequences[index])

	
print len(HAheaders), len(NAheaders)

i=0
while i<4:
	for i in range(0,len(HAsequencesDNA)):
		for j in range(0,len(HAsequencesDNA[i])-2,3):
			if HAsequencesDNA[i][j:j+3] == "TAA":
				HAsequencesDNA[i] = HAsequencesDNA[i][0:j]
			if HAsequencesDNA[i][i:i+3] == "TAG":
				HAsequencesDNA[i] = HAsequencesDNA[i][0:j]
			if HAsequencesDNA[i][i:i+3] == "TGA":
				HAsequencesDNA[i] = HAsequencesDNA[i][0:j]
		HAsequencesDNA[i] = HAsequencesDNA[i][3:]
			
	for i in range(0,len(NAsequencesDNA)):
		for j in range(0,len(NAsequencesDNA[i])-2,3):
			if NAsequencesDNA[i][j:j+3] == "TAA":
				NAsequencesDNA[i] = NAsequencesDNA[i][0:j]
			if NAsequencesDNA[i][i:i+3] == "TAG":
				NAsequencesDNA[i] = NAsequencesDNA[i][0:j]
			if NAsequencesDNA[i][i:i+3] == "TGA":
				NAsequencesDNA[i] = NAsequencesDNA[i][0:j]
		NAsequencesDNA[i] = NAsequencesDNA[i][3:]
		
	i+=1 
	
	
i=0
while i<4:
	for sequence in HAsequencesDNA:
		if sorted(set(sequence)) != ['A', 'C', 'G', 'T']:
			index = HAsequencesDNA.index(sequence)
			del HAheaders[index]
			del HAsequencesDNA[index]
			del NAheaders[index]
			del NAsequencesDNA[index]
			
	for sequence in NAsequencesDNA:
		if sorted(set(sequence)) != ['A', 'C', 'G', 'T']:
			index = NAsequencesDNA.index(sequence)
			del HAheaders[index]
			del HAsequencesDNA[index]
			del NAheaders[index]
			del NAsequencesDNA[index]
	
	i+=1
	
i=0
while i<4:
	for sequence in HAsequencesDNA:
		if len(sequence)<100:
			index = HAsequencesDNA.index(sequence)
			del HAheaders[index]
			del HAsequencesDNA[index]
			del NAheaders[index]
			del NAsequencesDNA[index]
			
	for sequence in NAsequencesDNA:
		if len(sequence)<100:
			index = NAsequencesDNA.index(sequence)
			del HAheaders[index]
			del HAsequencesDNA[index]
			del NAheaders[index]
			del NAsequencesDNA[index]
	
	i+=1

	
HAanalysisDNA = []
NAanalysisDNA = []

for sequence in HAsequencesDNA:
	HAanalysisDNA.append(DNA.dnaProfile(sequence))
	
for sequence in NAsequencesDNA:
	NAanalysisDNA.append(DNA.dnaProfile(sequence))
	
	
CpGDataHA = [] 
GCcontentHA = []
for y in HAanalysisDNA:
	GCcontentHA.append(y[0])
	CpGDataHA.append(y[1])
	
CpGDataNA = [] 
GCcontentNA = []
for y in NAanalysisDNA:
	GCcontentNA.append(y[0])
	CpGDataNA.append(y[1])


codons = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L","TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S","TAT":"Y", "TAC":"Y", "TAA":"", "TAG":"", "TGT":"C", "TGC":"C", "TGA":"", "TGG":"W","CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L","CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q","CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R","ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M","ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T","AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K","AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R","GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V","GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A","GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E","GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

HAsequencesProtein = []
for HA in HAsequencesDNA:
	psequence = ""
	for i in range(0,len(HA)-2,3):
		psequence += codons[HA[i:i+3]]
	HAsequencesProtein.append(psequence)
	
	
aminoacid_frequenciesHA = []
dipeptide_frequenciesHA = []
GRAVY_scoresHA = []
	
for sequence in HAsequencesProtein:
	amino_acids = {"G":0, "A":0, "L":0, "M":0, "F":0, "W":0, "K":0, "Q":0, "E":0, "S":0, "P":0, "V":0, "I":0, "C":0, "Y":0, "H":0, "R":0, "N":0, "D":0, "T":0}
	aminoacid_frequency = []

	for i in range(0,len(sequence)):
		amino_acids[sequence[i]] +=1.0

	i = 1
	for x,y in sorted(amino_acids.items()):
		aminoacid_frequency.append((i,(y/(len(sequence)))*100))
		i +=1
	
	aminoacid_frequenciesHA.append(aminoacid_frequency)
	
	dipeptides_dict = {'GW': 0, 'GV': 0, 'GT': 0, 'GS': 0, 'GR': 0, 'GQ': 0, 'GP': 0, 'GY': 0, 'GG': 0, 'GF': 0, 'GE': 0, 'GD': 0, 'GC': 0, 'GA': 0, 'GN': 0, 'GM': 0, 'GL': 0, 'GK': 0, 'GI': 0, 'GH': 0, 'ME': 0, 'MD': 0, 'MG': 0, 'MF': 0, 'MA': 0, 'MC': 0, 'MM': 0, 'ML': 0, 'MN': 0, 'MI': 0, 'MH': 0, 'MK': 0, 'MT': 0, 'MW': 0, 'MV': 0, 'MQ': 0, 'MP': 0, 'MS': 0, 'MR': 0, 'MY': 0, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'FA': 0, 'FC': 0, 'FD': 0, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0, 'SY': 0, 'SS': 0, 'SR': 0, 'SQ': 0, 'SP': 0, 'SW': 0, 'SV': 0, 'ST': 0, 'SK': 0, 'SI': 0, 'SH': 0, 'SN': 0, 'SM': 0, 'SL': 0, 'SC': 0, 'SA': 0, 'SG': 0, 'SF': 0, 'SE': 0, 'SD': 0, 'YI': 0, 'YH': 0, 'YK': 0, 'YM': 0, 'YL': 0, 'YN': 0, 'YA': 0, 'YC': 0, 'YE': 0, 'YD': 0, 'YG': 0, 'YF': 0, 'YY': 0, 'YQ': 0, 'YP': 0, 'YS': 0, 'YR': 0, 'YT': 0, 'YW': 0, 'YV': 0, 'LF': 0, 'LG': 0, 'LD': 0, 'LE': 0, 'LC': 0, 'LA': 0, 'LN': 0, 'LL': 0, 'LM': 0, 'LK': 0, 'LH': 0, 'LI': 0, 'LV': 0, 'LW': 0, 'LT': 0, 'LR': 0, 'LS': 0, 'LP': 0, 'LQ': 0, 'LY': 0, 'RT': 0, 'RV': 0, 'RW': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0, 'RY': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RA': 0, 'RC': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'VH': 0, 'IP': 0, 'EM': 0, 'EL': 0, 'EN': 0, 'EI': 0, 'EH': 0, 'EK': 0, 'EE': 0, 'ED': 0, 'EG': 0, 'EF': 0, 'EA': 0, 'EC': 0, 'IT': 0, 'EY': 0, 'VN': 0, 'ET': 0, 'EW': 0, 'EV': 0, 'EQ': 0, 'EP': 0, 'ES': 0, 'ER': 0, 'II': 0, 'VQ': 0, 'IK': 0, 'VT': 0, 'IN': 0, 'KC': 0, 'KA': 0, 'KG': 0, 'KF': 0, 'KE': 0, 'KD': 0, 'KK': 0, 'KI': 0, 'KH': 0, 'KN': 0, 'KM': 0, 'KL': 0, 'KS': 0, 'KR': 0, 'KQ': 0, 'KP': 0, 'KW': 0, 'KV': 0, 'KT': 0, 'KY': 0, 'DN': 0, 'DL': 0, 'DM': 0, 'DK': 0, 'DH': 0, 'DI': 0, 'DF': 0, 'DG': 0, 'DD': 0, 'DE': 0, 'DC': 0, 'DA': 0, 'DY': 0, 'DV': 0, 'DW': 0, 'DT': 0, 'DR': 0, 'DS': 0, 'DP': 0, 'DQ': 0, 'QQ': 0, 'QP': 0, 'QS': 0, 'QR': 0, 'QT': 0, 'QW': 0, 'QV': 0, 'QY': 0, 'QA': 0, 'QC': 0, 'QE': 0, 'QD': 0, 'QG': 0, 'QF': 0, 'QI': 0, 'QH': 0, 'QK': 0, 'QM': 0, 'QL': 0, 'QN': 0, 'WG': 0, 'WF': 0, 'WE': 0, 'WD': 0, 'WC': 0, 'WA': 0, 'WN': 0, 'WM': 0, 'WL': 0, 'WK': 0, 'WI': 0, 'WH': 0, 'WW': 0, 'WV': 0, 'WT': 0, 'WS': 0, 'WR': 0, 'WQ': 0, 'WP': 0, 'WY': 0, 'PR': 0, 'PS': 0, 'PP': 0, 'PQ': 0, 'PV': 0, 'PW': 0, 'PT': 0, 'PY': 0, 'PC': 0, 'PA': 0, 'PF': 0, 'PG': 0, 'PD': 0, 'PE': 0, 'PK': 0, 'PH': 0, 'PI': 0, 'PN': 0, 'PL': 0, 'PM': 0, 'CK': 0, 'CI': 0, 'CH': 0, 'CN': 0, 'CM': 0, 'CL': 0, 'CC': 0, 'CA': 0, 'CG': 0, 'CF': 0, 'CE': 0, 'CD': 0, 'CY': 0, 'CS': 0, 'CR': 0, 'CQ': 0, 'CP': 0, 'CW': 0, 'CV': 0, 'CT': 0, 'IY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0, 'VF': 0, 'VG': 0, 'IQ': 0, 'VI': 0, 'IS': 0, 'IR': 0, 'VL': 0, 'VM': 0, 'IW': 0, 'IV': 0, 'VP': 0, 'IH': 0, 'VR': 0, 'VS': 0, 'IM': 0, 'IL': 0, 'VV': 0, 'VW': 0, 'IA': 0, 'VY': 0, 'IC': 0, 'IE': 0, 'ID': 0, 'IG': 0, 'IF': 0, 'HY': 0, 'HR': 0, 'HS': 0, 'HP': 0, 'HQ': 0, 'HV': 0, 'HW': 0, 'HT': 0, 'HK': 0, 'HH': 0, 'HI': 0, 'HN': 0, 'HL': 0, 'HM': 0, 'HC': 0, 'HA': 0, 'HF': 0, 'HG': 0, 'HD': 0, 'HE': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0, 'NN': 0, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0, 'NG': 0, 'NY': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'TY': 0, 'TV': 0, 'TW': 0, 'TT': 0, 'TR': 0, 'TS': 0, 'TP': 0, 'TQ': 0, 'TN': 0, 'TL': 0, 'TM': 0, 'TK': 0, 'TH': 0, 'TI': 0, 'TF': 0, 'TG': 0, 'TD': 0, 'TE': 0, 'TC': 0, 'TA': 0, 'AA': 0, 'AC': 0, 'AE': 0, 'AD': 0, 'AG': 0, 'AF': 0, 'AI': 0, 'AH': 0, 'AK': 0, 'AM': 0, 'AL': 0, 'AN': 0, 'AQ': 0, 'AP': 0, 'AS': 0, 'AR': 0, 'AT': 0, 'AW': 0, 'AV': 0, 'AY': 0, 'VK': 0}
	dipeptide_frequency = []

	for i in range(0,len(sequence)-1):
		dipeptides_dict[sequence[i:i+2]] +=1.0

	strain_df = []
	i = 1
	for x,y in sorted(dipeptides_dict.items()):
		strain_df.append((i,(y/(len(sequence)-1))*100))
		i += 1

	dipeptide_frequenciesHA.append(strain_df)	
	
	
	hydropathy = {"A":1.800,"R":-4.500,"N":-3.500,"D":-3.500,"C":2.500,"Q":-3.500,"E":-3.500,"G":-0.400,"H":-3.200,"I":4.500,"L":3.800,"K":-3.900,"M":1.900,"F":2.800,"P":-1.600,"S":-0.800,"T":-0.700,"W":-0.900,"Y":-1.300,"V":4.200}
	
	GRAVYscore = 0
	for acid in sequence:
		GRAVYscore += hydropathy[acid]
	
	GRAVY_scoresHA.append(GRAVYscore)
	
print len(amino_acids)
print len(dipeptides_dict)	

NAsequencesProtein = []
for NA in NAsequencesDNA:
	psequence = ""
	for i in range(0,len(NA)-2,3):
		psequence += codons[NA[i:i+3]]
	NAsequencesProtein.append(psequence)
	
aminoacid_frequenciesNA = []
dipeptide_frequenciesNA = []
GRAVY_scoresNA = []
	
for sequence in NAsequencesProtein:
	amino_acids = {"G":0, "A":0, "L":0, "M":0, "F":0, "W":0, "K":0, "Q":0, "E":0, "S":0, "P":0, "V":0, "I":0, "C":0, "Y":0, "H":0, "R":0, "N":0, "D":0, "T":0}
	aminoacid_frequency = []

	for i in range(0,len(sequence)):
		amino_acids[sequence[i]] +=1.0

	i = 1
	for x,y in sorted(amino_acids.items()):
		aminoacid_frequency.append((i,(y/(len(sequence)))*100))
		i +=1
	
	aminoacid_frequenciesNA.append(aminoacid_frequency)
	
	dipeptides_dict = {'GW': 0, 'GV': 0, 'GT': 0, 'GS': 0, 'GR': 0, 'GQ': 0, 'GP': 0, 'GY': 0, 'GG': 0, 'GF': 0, 'GE': 0, 'GD': 0, 'GC': 0, 'GA': 0, 'GN': 0, 'GM': 0, 'GL': 0, 'GK': 0, 'GI': 0, 'GH': 0, 'ME': 0, 'MD': 0, 'MG': 0, 'MF': 0, 'MA': 0, 'MC': 0, 'MM': 0, 'ML': 0, 'MN': 0, 'MI': 0, 'MH': 0, 'MK': 0, 'MT': 0, 'MW': 0, 'MV': 0, 'MQ': 0, 'MP': 0, 'MS': 0, 'MR': 0, 'MY': 0, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'FA': 0, 'FC': 0, 'FD': 0, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0, 'SY': 0, 'SS': 0, 'SR': 0, 'SQ': 0, 'SP': 0, 'SW': 0, 'SV': 0, 'ST': 0, 'SK': 0, 'SI': 0, 'SH': 0, 'SN': 0, 'SM': 0, 'SL': 0, 'SC': 0, 'SA': 0, 'SG': 0, 'SF': 0, 'SE': 0, 'SD': 0, 'YI': 0, 'YH': 0, 'YK': 0, 'YM': 0, 'YL': 0, 'YN': 0, 'YA': 0, 'YC': 0, 'YE': 0, 'YD': 0, 'YG': 0, 'YF': 0, 'YY': 0, 'YQ': 0, 'YP': 0, 'YS': 0, 'YR': 0, 'YT': 0, 'YW': 0, 'YV': 0, 'LF': 0, 'LG': 0, 'LD': 0, 'LE': 0, 'LC': 0, 'LA': 0, 'LN': 0, 'LL': 0, 'LM': 0, 'LK': 0, 'LH': 0, 'LI': 0, 'LV': 0, 'LW': 0, 'LT': 0, 'LR': 0, 'LS': 0, 'LP': 0, 'LQ': 0, 'LY': 0, 'RT': 0, 'RV': 0, 'RW': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0, 'RY': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RA': 0, 'RC': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'VH': 0, 'IP': 0, 'EM': 0, 'EL': 0, 'EN': 0, 'EI': 0, 'EH': 0, 'EK': 0, 'EE': 0, 'ED': 0, 'EG': 0, 'EF': 0, 'EA': 0, 'EC': 0, 'IT': 0, 'EY': 0, 'VN': 0, 'ET': 0, 'EW': 0, 'EV': 0, 'EQ': 0, 'EP': 0, 'ES': 0, 'ER': 0, 'II': 0, 'VQ': 0, 'IK': 0, 'VT': 0, 'IN': 0, 'KC': 0, 'KA': 0, 'KG': 0, 'KF': 0, 'KE': 0, 'KD': 0, 'KK': 0, 'KI': 0, 'KH': 0, 'KN': 0, 'KM': 0, 'KL': 0, 'KS': 0, 'KR': 0, 'KQ': 0, 'KP': 0, 'KW': 0, 'KV': 0, 'KT': 0, 'KY': 0, 'DN': 0, 'DL': 0, 'DM': 0, 'DK': 0, 'DH': 0, 'DI': 0, 'DF': 0, 'DG': 0, 'DD': 0, 'DE': 0, 'DC': 0, 'DA': 0, 'DY': 0, 'DV': 0, 'DW': 0, 'DT': 0, 'DR': 0, 'DS': 0, 'DP': 0, 'DQ': 0, 'QQ': 0, 'QP': 0, 'QS': 0, 'QR': 0, 'QT': 0, 'QW': 0, 'QV': 0, 'QY': 0, 'QA': 0, 'QC': 0, 'QE': 0, 'QD': 0, 'QG': 0, 'QF': 0, 'QI': 0, 'QH': 0, 'QK': 0, 'QM': 0, 'QL': 0, 'QN': 0, 'WG': 0, 'WF': 0, 'WE': 0, 'WD': 0, 'WC': 0, 'WA': 0, 'WN': 0, 'WM': 0, 'WL': 0, 'WK': 0, 'WI': 0, 'WH': 0, 'WW': 0, 'WV': 0, 'WT': 0, 'WS': 0, 'WR': 0, 'WQ': 0, 'WP': 0, 'WY': 0, 'PR': 0, 'PS': 0, 'PP': 0, 'PQ': 0, 'PV': 0, 'PW': 0, 'PT': 0, 'PY': 0, 'PC': 0, 'PA': 0, 'PF': 0, 'PG': 0, 'PD': 0, 'PE': 0, 'PK': 0, 'PH': 0, 'PI': 0, 'PN': 0, 'PL': 0, 'PM': 0, 'CK': 0, 'CI': 0, 'CH': 0, 'CN': 0, 'CM': 0, 'CL': 0, 'CC': 0, 'CA': 0, 'CG': 0, 'CF': 0, 'CE': 0, 'CD': 0, 'CY': 0, 'CS': 0, 'CR': 0, 'CQ': 0, 'CP': 0, 'CW': 0, 'CV': 0, 'CT': 0, 'IY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0, 'VF': 0, 'VG': 0, 'IQ': 0, 'VI': 0, 'IS': 0, 'IR': 0, 'VL': 0, 'VM': 0, 'IW': 0, 'IV': 0, 'VP': 0, 'IH': 0, 'VR': 0, 'VS': 0, 'IM': 0, 'IL': 0, 'VV': 0, 'VW': 0, 'IA': 0, 'VY': 0, 'IC': 0, 'IE': 0, 'ID': 0, 'IG': 0, 'IF': 0, 'HY': 0, 'HR': 0, 'HS': 0, 'HP': 0, 'HQ': 0, 'HV': 0, 'HW': 0, 'HT': 0, 'HK': 0, 'HH': 0, 'HI': 0, 'HN': 0, 'HL': 0, 'HM': 0, 'HC': 0, 'HA': 0, 'HF': 0, 'HG': 0, 'HD': 0, 'HE': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0, 'NN': 0, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0, 'NG': 0, 'NY': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'TY': 0, 'TV': 0, 'TW': 0, 'TT': 0, 'TR': 0, 'TS': 0, 'TP': 0, 'TQ': 0, 'TN': 0, 'TL': 0, 'TM': 0, 'TK': 0, 'TH': 0, 'TI': 0, 'TF': 0, 'TG': 0, 'TD': 0, 'TE': 0, 'TC': 0, 'TA': 0, 'AA': 0, 'AC': 0, 'AE': 0, 'AD': 0, 'AG': 0, 'AF': 0, 'AI': 0, 'AH': 0, 'AK': 0, 'AM': 0, 'AL': 0, 'AN': 0, 'AQ': 0, 'AP': 0, 'AS': 0, 'AR': 0, 'AT': 0, 'AW': 0, 'AV': 0, 'AY': 0, 'VK': 0}
	dipeptide_frequency = []

	for i in range(0,len(sequence)-1):
		dipeptides_dict[sequence[i:i+2]] +=1.0

	strain_df = []
	i = 1
	for x,y in sorted(dipeptides_dict.items()):
		strain_df.append((i,(y/(len(sequence)-1))*100))
		i += 1

	dipeptide_frequenciesNA.append(strain_df)	
	
	
	hydropathy = {"A":1.800,"R":-4.500,"N":-3.500,"D":-3.500,"C":2.500,"Q":-3.500,"E":-3.500,"G":-0.400,"H":-3.200,"I":4.500,"L":3.800,"K":-3.900,"M":1.900,"F":2.800,"P":-1.600,"S":-0.800,"T":-0.700,"W":-0.900,"Y":-1.300,"V":4.200}
	
	GRAVYscore = 0
	for acid in sequence:
		GRAVYscore += hydropathy[acid]
	
	GRAVY_scoresNA.append(GRAVYscore)
	
FinalLists = []

for i in range(0,len(HAheaders)):
	signatures = []
	signatures.append(GCcontentHA[i])
	signatures.append(CpGDataHA[i])
	signatures.append(GRAVY_scoresHA[i])
	for x,y in sorted(aminoacid_frequenciesHA[i]):
		signatures.append(y)
	for x,y in sorted(dipeptide_frequenciesHA[i]):
		signatures.append(y)
	signatures.append(GCcontentNA[i])
	signatures.append(CpGDataNA[i])
	signatures.append(GRAVY_scoresNA[i])
	for x,y in sorted(aminoacid_frequenciesNA[i]):
		signatures.append(y)
	for x,y in sorted(dipeptide_frequenciesNA[i]):
		signatures.append(y)
	FinalLists.append(signatures)
	
AverageStrain = []
for j in range(0,len(FinalLists[0])):
	sum = 0
	for i in range(0,len(FinalLists)):
		sum += FinalLists[i][j]
	average = sum/len(FinalLists)
	AverageStrain.append(average)
	
print AverageStrain 

#top 5
distances = []
for i in range(0,len(FinalLists)):
	sum = 0
	for j in range(0,len(FinalLists[i])):
		sum += ((FinalLists[i][j]) - AverageStrain[j])**2
	distance = math.sqrt(sum)
	distances.append(distance)

closest_neighbor = distances.index(min(distances))	
print HAheaders[closest_neighbor]
print NAheaders[closest_neighbor]

print HAsequencesProtein[closest_neighbor]
	
	
outputCSV = open(outputAnalysis, "w")
outputCSV.writelines(",GC content, CpG ratio, GRAVY Score,")
for acid in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
	outputCSV.writelines(acid + ",")
for dipeptide in ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL', 'AM', 'AN', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AV', 'AW', 'AY', 'CA', 'CC', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT', 'CV', 'CW', 'CY', 'DA', 'DC', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN', 'DP', 'DQ', 'DR', 'DS', 'DT', 'DV', 'DW', 'DY', 'EA', 'EC', 'ED', 'EE', 'EF', 'EG', 'EH', 'EI', 'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER', 'ES', 'ET', 'EV', 'EW', 'EY', 'FA', 'FC', 'FD', 'FE', 'FF', 'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT', 'FV', 'FW', 'FY', 'GA', 'GC', 'GD', 'GE', 'GF', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS', 'GT', 'GV', 'GW', 'GY', 'HA', 'HC', 'HD', 'HE', 'HF', 'HG', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP', 'HQ', 'HR', 'HS', 'HT', 'HV', 'HW', 'HY', 'IA', 'IC', 'ID', 'IE', 'IF', 'IG', 'IH', 'II', 'IK', 'IL', 'IM', 'IN', 'IP', 'IQ', 'IR', 'IS', 'IT', 'IV', 'IW', 'IY', 'KA', 'KC', 'KD', 'KE', 'KF', 'KG', 'KH', 'KI', 'KK', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR', 'KS', 'KT', 'KV', 'KW', 'KY', 'LA', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LK', 'LL', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT', 'LV', 'LW', 'LY', 'MA', 'MC', 'MD', 'ME', 'MF', 'MG', 'MH', 'MI', 'MK', 'ML', 'MM', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY', 'NA', 'NC', 'ND', 'NE', 'NF', 'NG', 'NH', 'NI', 'NK', 'NL', 'NM', 'NN', 'NP', 'NQ', 'NR', 'NS', 'NT', 'NV', 'NW', 'NY', 'PA', 'PC', 'PD', 'PE', 'PF', 'PG', 'PH', 'PI', 'PK', 'PL', 'PM', 'PN', 'PP', 'PQ', 'PR', 'PS', 'PT', 'PV', 'PW', 'PY', 'QA', 'QC', 'QD', 'QE', 'QF', 'QG', 'QH', 'QI', 'QK', 'QL', 'QM', 'QN', 'QP', 'QQ', 'QR', 'QS', 'QT', 'QV', 'QW', 'QY', 'RA', 'RC', 'RD', 'RE', 'RF', 'RG', 'RH', 'RI', 'RK', 'RL', 'RM', 'RN', 'RP', 'RQ', 'RR', 'RS', 'RT', 'RV', 'RW', 'RY', 'SA', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SK', 'SL', 'SM', 'SN', 'SP', 'SQ', 'SR', 'SS', 'ST', 'SV', 'SW', 'SY', 'TA', 'TC', 'TD', 'TE', 'TF', 'TG', 'TH', 'TI', 'TK', 'TL', 'TM', 'TN', 'TP', 'TQ', 'TR', 'TS', 'TT', 'TV', 'TW', 'TY', 'VA', 'VC', 'VD', 'VE', 'VF', 'VG', 'VH', 'VI', 'VK', 'VL', 'VM', 'VN', 'VP', 'VQ', 'VR', 'VS', 'VT', 'VV', 'VW', 'VY', 'WA', 'WC', 'WD', 'WE', 'WF', 'WG', 'WH', 'WI', 'WK', 'WL', 'WM', 'WN', 'WP', 'WQ', 'WR', 'WS', 'WT', 'WV', 'WW', 'WY', 'YA', 'YC', 'YD', 'YE', 'YF', 'YG', 'YH', 'YI', 'YK', 'YL', 'YM', 'YN', 'YP', 'YQ', 'YR', 'YS', 'YT', 'YV', 'YW', 'YY']:
	outputCSV.writelines(dipeptide + ",")
outputCSV.writelines("GC content, CpG ratio, GRAVY Score,")
for acid in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
	outputCSV.writelines(acid + ",")
for dipeptide in ['AA', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AK', 'AL', 'AM', 'AN', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AV', 'AW', 'AY', 'CA', 'CC', 'CD', 'CE', 'CF', 'CG', 'CH', 'CI', 'CK', 'CL', 'CM', 'CN', 'CP', 'CQ', 'CR', 'CS', 'CT', 'CV', 'CW', 'CY', 'DA', 'DC', 'DD', 'DE', 'DF', 'DG', 'DH', 'DI', 'DK', 'DL', 'DM', 'DN', 'DP', 'DQ', 'DR', 'DS', 'DT', 'DV', 'DW', 'DY', 'EA', 'EC', 'ED', 'EE', 'EF', 'EG', 'EH', 'EI', 'EK', 'EL', 'EM', 'EN', 'EP', 'EQ', 'ER', 'ES', 'ET', 'EV', 'EW', 'EY', 'FA', 'FC', 'FD', 'FE', 'FF', 'FG', 'FH', 'FI', 'FK', 'FL', 'FM', 'FN', 'FP', 'FQ', 'FR', 'FS', 'FT', 'FV', 'FW', 'FY', 'GA', 'GC', 'GD', 'GE', 'GF', 'GG', 'GH', 'GI', 'GK', 'GL', 'GM', 'GN', 'GP', 'GQ', 'GR', 'GS', 'GT', 'GV', 'GW', 'GY', 'HA', 'HC', 'HD', 'HE', 'HF', 'HG', 'HH', 'HI', 'HK', 'HL', 'HM', 'HN', 'HP', 'HQ', 'HR', 'HS', 'HT', 'HV', 'HW', 'HY', 'IA', 'IC', 'ID', 'IE', 'IF', 'IG', 'IH', 'II', 'IK', 'IL', 'IM', 'IN', 'IP', 'IQ', 'IR', 'IS', 'IT', 'IV', 'IW', 'IY', 'KA', 'KC', 'KD', 'KE', 'KF', 'KG', 'KH', 'KI', 'KK', 'KL', 'KM', 'KN', 'KP', 'KQ', 'KR', 'KS', 'KT', 'KV', 'KW', 'KY', 'LA', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LK', 'LL', 'LM', 'LN', 'LP', 'LQ', 'LR', 'LS', 'LT', 'LV', 'LW', 'LY', 'MA', 'MC', 'MD', 'ME', 'MF', 'MG', 'MH', 'MI', 'MK', 'ML', 'MM', 'MN', 'MP', 'MQ', 'MR', 'MS', 'MT', 'MV', 'MW', 'MY', 'NA', 'NC', 'ND', 'NE', 'NF', 'NG', 'NH', 'NI', 'NK', 'NL', 'NM', 'NN', 'NP', 'NQ', 'NR', 'NS', 'NT', 'NV', 'NW', 'NY', 'PA', 'PC', 'PD', 'PE', 'PF', 'PG', 'PH', 'PI', 'PK', 'PL', 'PM', 'PN', 'PP', 'PQ', 'PR', 'PS', 'PT', 'PV', 'PW', 'PY', 'QA', 'QC', 'QD', 'QE', 'QF', 'QG', 'QH', 'QI', 'QK', 'QL', 'QM', 'QN', 'QP', 'QQ', 'QR', 'QS', 'QT', 'QV', 'QW', 'QY', 'RA', 'RC', 'RD', 'RE', 'RF', 'RG', 'RH', 'RI', 'RK', 'RL', 'RM', 'RN', 'RP', 'RQ', 'RR', 'RS', 'RT', 'RV', 'RW', 'RY', 'SA', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH', 'SI', 'SK', 'SL', 'SM', 'SN', 'SP', 'SQ', 'SR', 'SS', 'ST', 'SV', 'SW', 'SY', 'TA', 'TC', 'TD', 'TE', 'TF', 'TG', 'TH', 'TI', 'TK', 'TL', 'TM', 'TN', 'TP', 'TQ', 'TR', 'TS', 'TT', 'TV', 'TW', 'TY', 'VA', 'VC', 'VD', 'VE', 'VF', 'VG', 'VH', 'VI', 'VK', 'VL', 'VM', 'VN', 'VP', 'VQ', 'VR', 'VS', 'VT', 'VV', 'VW', 'VY', 'WA', 'WC', 'WD', 'WE', 'WF', 'WG', 'WH', 'WI', 'WK', 'WL', 'WM', 'WN', 'WP', 'WQ', 'WR', 'WS', 'WT', 'WV', 'WW', 'WY', 'YA', 'YC', 'YD', 'YE', 'YF', 'YG', 'YH', 'YI', 'YK', 'YL', 'YM', 'YN', 'YP', 'YQ', 'YR', 'YS', 'YT', 'YV', 'YW', 'YY']:
	outputCSV.writelines(dipeptide + ",")
outputCSV.writelines("\n")

for i in range(0,len(GCcontentHA)):
	outputCSV.writelines(HAheaders[i] + ",")
	for value in FinalLists[i]:
		outputCSV.writelines(str(value) + ",")
	outputCSV.writelines("\n")
	
outputCSV.writelines(HAsequencesProtein[closest_neighbor])
outputCSV.writelines("\n")
outputCSV.writelines(NAsequencesProtein[closest_neighbor])
outputCSV.writelines("\n")
outputCSV.writelines(HAheaders[closest_neighbor])
	
outputCSV.close()
