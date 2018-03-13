import DNA_calculators as dna	
import sklearn
from sklearn.cluster import KMeans
import matplotlib.pyplot as plots
import numpy
import time
import sys

inputSequences = sys.argv[1] 
outputAnalysis = sys.argv[2] 

#opens and reads in the input file
fsequences = open(inputSequences, "r")
sequencesData = fsequences.readlines()
fsequences.close()

headers = []
sequences = []
gene = ""
for line in sequencesData:
	if line[0] == ">":
		headers.append(line)
		if len(gene) != 0:
			sequences.append(gene)
			gene = ""
	else:
		gene += line.strip()		
sequences.append(gene)

i=0
while i<3:
	for sequence in sequences:
		if sorted(set(sequence)) != ['A', 'C', 'G', 'T']:
			index = sequences.index(sequence)
			del sequences[index]
			del headers[index]
	i+=1	

strainProfiles = {}
for i in range(0,len(sequences)):
	strainProfiles[headers[i]] = dna.dnaProfile(sequences[i])
	
CpGData = [] 
GCratio = []
AT_GCratio = []
dinucleotide = []
codon_usage = []
gene_numbers = []
for i in range(0,len(strainProfiles)):
	gene_numbers.append(i)

for x,y in strainProfiles.items():
	CpGData.append(y[2])
	GCratio.append(y[0])
	AT_GCratio.append(y[1])
	dinucleotide.append(y[3])
	codon_usage.append(y[4])

geneNumber = [range(0,len(CpGData))]
	
plots.scatter(GCratio, CpGData)
plots.show()

CpGDatapoints = []
for i in range(0,len(CpGData)):
	CpGDatapoints.append((GCratio[i], CpGData[i]))
fluClusters = KMeans(n_clusters = 3).fit(CpGDatapoints)
print fluClusters.cluster_centers_



outputInfo = open(outputAnalysis, "w")

outputInfo.writelines("GC content, AT/GC ratio, CpG ratio, dinucleotide frequency, codon frequency")
outputInfo.writelines("\n\n")
for x,y in strainProfiles.items():
		outputInfo.writelines(x)
		outputInfo.writelines(str(y))
		outputInfo.writelines("\n\n")

outputInfo.close()
