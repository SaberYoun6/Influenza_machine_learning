import sys




class file_reader(object):
	def __init__(files):
		self.files =files
	def file_readers():
		fsequences = open(self.files,'r')
		sequencesData = fsequences.readlines()
		fsequences.close()
		return sequencesData
	def header_setter(sequencesD):
		headers = []
		for line in sequencesD:
			if ">" in line:
				headers.append(line)
		return headers
	def seq_setter(sequencesD):
		sequences = []
		gene = ""
		for i in sequencesD:
			if len(gene) != 0:
				sequences.append(gene)
				gene = ""
			else:
				gene += line.strip()		
			sequences.append(gene)
		return sequences
	def NA_sequence(header,sequences):
		NAheaders=[]
		NAsequences = []
		for header in headers:
			if "neuraminidase" in header:
				NAheaders.append(header)
				index = headers.index(header)
				NAsequences.append(sequences[index])
		return (NAheaders,NAsequences)

	def HA_sequence(header,sequences):
		HAheaders = []
		HAsequences = []
		for header in headers:
			if "hemagglutinin" in header:
				HAheaders.append(header)
				index = headers.index(header)
				HAsequences.append(sequences[index])
		return (HAheaders,HAsequences)


