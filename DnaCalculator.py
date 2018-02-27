import DNA_calculators as dna	
import numpy

#opens and reads in the input file
fsequences = open("flu.fasta", "r")
fluData = fsequences.readlines()
fsequences.close()

headers = []
sequences = []
gene = ""
for line in fluData:
	if line[0] == ">":
		headers.append(line)
		if len(gene) != 0:
			sequences.append(gene)
			gene = ""
	else:
		gene += line.strip()
		
sequences.append(gene)

strainProfiles = {}
for i in range(0,len(sequences)):
	strainProfiles[headers[i]] = dna.fullProfile(sequences[i])

outputInfo = open("preliminaryOutput.txt", "w")

outputInfo.writelines("CpG ratio, GC contenct, AT/GC ratio, dinucleotide frequency, codon frequency")
outputInfo.writelines("\n\n")
for x,y in strainProfiles.items():
		outputInfo.writelines(x)
		outputInfo.writelines(str(y))
		outputInfo.writelines("\n\n")

outputInfo.close()
