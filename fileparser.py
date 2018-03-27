fsequences = open("0910HANAH1N1Genes.fasta","r")
sequencesData = fsequences.readlines()
fsequences.close()


headers = []
sequences = []
gene = ""
for line in sequencesData:
	if ">" in line:
		headers.append(line)
		if len(gene) != 0:
			sequences.append(gene)
			gene = ""
	else:
		gene += line.strip()		
sequences.append(gene)

print len(headers)
print len(sequences)

HAheaders = []
HAsequences = []
NAheaders = []
NAsequences = []

for header in headers:
	if "neuraminidase" in header:
		NAheaders.append(header)
		index = headers.index(header)
		NAsequences.append(sequences[index])
	
	if "hemagglutinin" in header:
		HAheaders.append(header)
		index = headers.index(header)
		HAsequences.append(sequences[index])
		
print len(HAheaders)
print len(HAsequences)

print len(NAheaders)
print len(NAsequences)