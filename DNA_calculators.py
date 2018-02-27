#These are calculators for the DNA motifs that we'll be using. They will likely be defined functions with return values that will
#run on loops for each sequence in the influenza strains
#Because of the technique for collecting flu gene sequences, nucleotides
#have T instead of the U in flu RNA. 
#For order of nucleotides, we go in the order A,C,G,C
#the initial dictionaries are not in that order due to the mutability, but when 
#values from the dict are added to lists they are in the order A,C,G,T
#For example, AA, AC, AG, AT or AAA, AAC, AAG, AAT


#the fullProfile produces a list in the order CpG ratio, GC contenct, AT/GC ratio,
#dinucleotide frequency, and then codon frequency 
def fullProfile(sequence):
	length = float(len(sequence))
	strainProfile = []
	A, C, G, T = 0, 0, 0, 0
	for i in range(0,len(sequence)): 
		if sequence[i] == "A":
			A +=1
		if sequence[i] == "C":
			C +=1
		if sequence[i] == "T":
			T +=1
		if sequence[i] == "G":
			G +=1

	GC_content = (G+C)/length
	strainProfile.append(GC_content)

	AT_CG_ratio = float(A+T)/float(C+G)
	strainProfile.append(AT_CG_ratio)
	
	CpGcounter = 0.0
	for i in range(0,len(sequence)-1):
		if sequence[i] == "C" and i != len(sequence):
			if sequence[i+1] =="G":
				CpGcounter += 1
		else:
			continue

	Gfrequency = G/length
	Cfrequency = C/length
	CpGexpected = int((C*G)*length)

	CpGRatio = CpGcounter/CpGexpected
	strainProfile.append(CpGRatio)


	dinucleotide_dict = {'AA': 0, 'AC': 0, 'GT': 0, 'AG': 0, 'CC': 0, 'CA': 0, 'CG': 0, 'TT': 0, 'GG': 0, 'GC': 0, 'AT': 0, 'GA': 0, 'TG': 0, 'TA': 0, 'TC': 0, 'CT': 0}
	dinucleotide_frequency = []
	for i in range(0,len(sequence)-1):
		dinucleotide_dict[sequence[i:i+2]] += 1

	for x,y in sorted(dinucleotide_dict.items()):
		dinucleotide_frequency.append(y/length)
	
	strainProfile.append(dinucleotide_frequency)

	codons = {"ATT":0, "ATC":0, "ATA":0, "CTT":0, "CTC":0, "CTA":0, "CTG":0, "TTA":0, "TTG":0, "GTT":0, "GTC":0, "GTA":0, "GTG":0, "TTT":0, "TTC":0, "ATG":0, "TGT":0, "TGC":0, "GCT":0, "GCC":0, "GCA":0, "GCG":0, "GGT":0, "GGC":0, "GGA":0, "GGG":0, "CCT":0, "CCC":0, "CCA":0, "CCG":0, "ACT":0, "ACC":0, "ACA":0, "ACG":0, "TCT":0, "TCC":0, "TCA":0, "TCG":0, "AGT":0, "AGC":0, "TAT":0, "TAC":0, "TGG":0, "CAA":0, "CAG":0, "AAT":0, "AAC":0, "CAT":0, "CAC":0, "GAA":0, "GAG":0, "GAT":0, "GAC":0, "AAA":0, "AAG":0, "CGT":0, "CGC":0, "CGA":0, "CGG":0, "AGA":0, "AGG":0, "TAA":0, "TAG":0, "TGA":0} 
	codon_frequency = []
	for i in range(0, len(sequence)-2):
		codons[sequence[i:i+3]] +=1

	for x,y in sorted(codons.items()):
		codon_frequency.append(y/length)
	strainProfile.append(codon_frequency)
	
	return strainProfile