#These are calculators for the DNA motifs that we'll be using. They will likely be defined functions with return values that will
#run on loops for each sequence in the influenza strains
#For order of nucleotides, we go in the order A,C,T,G

sequence = "ATCGTCGATCGTAGCTGATCGTAGTTTGGGCCAAGCTGCTGCCGCGCGCTGATCGCTCGCTAGCTCGT"
CpG_ratios, dinucleotide_frequencies, codon_frequencies, GC_content, AT_CG_ratios = [], [], [], [], []
length = float(len(sequence))


A, C, T, G = 0, 0, 0, 0
for i in range(0,len(sequence)): 
	if sequence[i] == "A":
		A +=1
	if sequence[i] == "C":
		C +=1
	if sequence[i] == "T":
		T +=1
	if sequence[i] == "G":
		G +=1

GC_content.append((G+C)/length)
AT_CG_ratios.append((float(A+T))/C+G)
	
CpGcounter = 0
for i in range(0,len(sequence)-1):
	if sequence[i] == "C" and i != len(sequence):
		if sequence[i+1] =="G":
			CpGcounter += 1
	else:
		continue

Gfrequency = G/length
Cfrequency = C/length
CpGexpected = int((C*G)*length)

CpGRatio = CpG_counter/CpGexpected
CpG_ratios.append(CpGRatio)




AA, AC, AT, AG, CA, CC, CT, CG, TA, TC, TT, TG, GA, GC, GT, GG = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
for i in range(0,len(sequence)-1):
	
	
	if sequence[i] == "A":
		if sequence[i+1] == "A":
			AA +=1
		
		if sequence[i+1] == "C":
			AC +=1
		
		if sequence[i+1] == "T":
			AT +=1
		
		if sequence[i+1] == "G":
			AG +=1
	
	if sequence[i] == "C":
		if sequence[i+1] == "A":
			CA +=1
		
		if sequence[i+1] == "C":
			CC +=1
		
		if sequence[i+1] == "T":
			CT +=1
		
		if sequence[i+1] == "G":
			CG +=1
	
	if sequence[i] == "T":
		if sequence[i+1] == "A":
			TA +=1
		
		if sequence[i+1] == "C":
			TC +=1
		
		if sequence[i+1] == "T":
			TT +=1
		
		if sequence[i+1] == "G":
			TG +=1
	
	if sequence[i] == "G":
		if sequence[i+1] == "A":
			GA +=1
		
		if sequence[i+1] == "C":
			GC +=1
		
		if sequence[i+1] == "T":
			GT +=1
		
		if sequence[i+1] == "G":
			GG +=1
	
dinucleotide_frequencies.append([AA/length, AC/length, AT/length, AG/length, CA/length, CC/length, CT/length, CG/length, TA/length, TC/length, TT/length, TG/length, GA/length, GC/length, GT/length, GG/length])
			
codons = {"ATT":0, "ATC":0, "ATA":0, "CTT":0, "CTC":0, "CTA":0, "CTG":0, "TTA":0, "TTG":0, "GTT":0, "GTC":0, "GTA":0, "GTG":0, "TTT":0, "TTC":0, "ATG":0, "TGT":0, "TGC":0, "GCT":0, "GCC":0, "GCA":0, "GCG":0, "GGT":0, "GGC":0, "GGA":0, "GGG":0, "CCT":0, "CCC":0, "CCA":0, "CCG":0, "ACT":0, "ACC":0, "ACA":0, "ACG":0, "TCT":0, "TCC":0, "TCA":0, "TCG":0, "AGT":0, "AGC":0, "TAT":0, "TAC":0, "TGG":0, "CAA":0, "CAG":0, "AAT":0, "AAC":0, "CAT":0, "CAC":0, "GAA":0, "GAG":0, "GAT":0, "GAC":0, "AAA":0, "AAG":0, "CGT":0, "CGC":0, "CGA":0, "CGG":0, "AGA":0, "AGG":0, "TAA":0, "TAG":0, "TGA":0} 
codon_frequency = []
for i in range(0, len(sequence)-2):
	codons[sequence[i:i+3]] +=1

for x,y in codons.items():
	codon_frequency.append(y/length)

codon_frequencies.append(codon_frequency)
	

	
			
print CpG_ratios, dinucleotide_frequencies, codon_frequencies, GC_content, AT_CG_ratios
