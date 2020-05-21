'''

'''




class file_reader(object):
	def __init__(self,files):
		self.files =files

	'''
	
	'''
	def header_setter(self,sequencesD):
		headers = []
		for line in sequencesD:
			if ">" in line:
				headers.append(line)
		return headers
	'''

	'''

	def seq_setter(self,sequencesD):
		sequences = []
		gene = ""
		
		for line in sequencesD:
			if line.startswith(">"):
				continue
			line.rstrip().lstrip("\n")
			if len(gene) != 0:
				gene+= line
				gene = ""
		sequences.append(gene)
		return sequences
	'''

	'''

	def file_readers(self):
		fsequences = open(self.files,'r')
		sequencesData = fsequences.readlines()
		fsequences.close()
		return  sequencesData
	'''


	'''
	def NA_sequence(self,head,seq):
		NAheaders=[]
		NAsequences = []
		for header in headers:
			if "neuraminidase" in header:
				NAheaders.append(header)
				index = headers.index(header)
				NAsequences.append(sequences[index])
		return (NAheaders,NAsequences)

	'''
	this was any ideal to try and split up headers by the Hemagglutinin and allow for the sequeces to be founds

	'''
	def HA_sequence(self,head,seq):
		HAheaders = []
		HAsequences = []
		for header in headers:
			if "hemagglutinin" in header:
				HAheaders.append(header)
				index = headers.index(header)
				HAsequences.append(sequences[index])
		return (HAheaders,HAsequences)


