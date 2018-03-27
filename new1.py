import random
import math
import matplotlib.pyplot as plots
	
#opens and reads in the input file
fsequences = open("hemagglutininProteins.fasta", "r")
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
		for acid in set(sequence):
			if acid in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
				continue
			else:
				index = sequences.index(sequence)
				del sequences[index]
				del headers[index]
				break
	i+=1	
	
	for sequence in sequences:
		if sequence[0] != "M":
			index = sequences.index(sequence)
			del sequences[index]
			del headers[index]
			break
	
aminoacid_frequencies = []
dipeptide_frequencies = []
	
	
for sequence in sequences:

	amino_acids = {"G":0, "A":0, "L":0, "M":0, "F":0, "W":0, "K":0, "Q":0, "E":0, "S":0, "P":0, "V":0, "I":0, "C":0, "Y":0, "H":0, "R":0, "N":0, "D":0, "T":0}
	aminoacid_frequency = []

	for i in range(0,len(sequence)):
		amino_acids[sequence[i]] +=1.0

	i = 1
	for x,y in sorted(amino_acids.items()):
		aminoacid_frequency.append((i,y))
		i +=1.0
	
	aminoacid_frequencies.append(aminoacid_frequency)
	
	dipeptides_dict = {'GW': 0, 'GV': 0, 'GT': 0, 'GS': 0, 'GR': 0, 'GQ': 0, 'GP': 0, 'GY': 0, 'GG': 0, 'GF': 0, 'GE': 0, 'GD': 0, 'GC': 0, 'GA': 0, 'GN': 0, 'GM': 0, 'GL': 0, 'GK': 0, 'GI': 0, 'GH': 0, 'ME': 0, 'MD': 0, 'MG': 0, 'MF': 0, 'MA': 0, 'MC': 0, 'MM': 0, 'ML': 0, 'MN': 0, 'MI': 0, 'MH': 0, 'MK': 0, 'MT': 0, 'MW': 0, 'MV': 0, 'MQ': 0, 'MP': 0, 'MS': 0, 'MR': 0, 'MY': 0, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'FA': 0, 'FC': 0, 'FD': 0, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0, 'SY': 0, 'SS': 0, 'SR': 0, 'SQ': 0, 'SP': 0, 'SW': 0, 'SV': 0, 'ST': 0, 'SK': 0, 'SI': 0, 'SH': 0, 'SN': 0, 'SM': 0, 'SL': 0, 'SC': 0, 'SA': 0, 'SG': 0, 'SF': 0, 'SE': 0, 'SD': 0, 'YI': 0, 'YH': 0, 'YK': 0, 'YM': 0, 'YL': 0, 'YN': 0, 'YA': 0, 'YC': 0, 'YE': 0, 'YD': 0, 'YG': 0, 'YF': 0, 'YY': 0, 'YQ': 0, 'YP': 0, 'YS': 0, 'YR': 0, 'YT': 0, 'YW': 0, 'YV': 0, 'LF': 0, 'LG': 0, 'LD': 0, 'LE': 0, 'LC': 0, 'LA': 0, 'LN': 0, 'LL': 0, 'LM': 0, 'LK': 0, 'LH': 0, 'LI': 0, 'LV': 0, 'LW': 0, 'LT': 0, 'LR': 0, 'LS': 0, 'LP': 0, 'LQ': 0, 'LY': 0, 'RT': 0, 'RV': 0, 'RW': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0, 'RY': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RA': 0, 'RC': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'VH': 0, 'IP': 0, 'EM': 0, 'EL': 0, 'EN': 0, 'EI': 0, 'EH': 0, 'EK': 0, 'EE': 0, 'ED': 0, 'EG': 0, 'EF': 0, 'EA': 0, 'EC': 0, 'IT': 0, 'EY': 0, 'VN': 0, 'ET': 0, 'EW': 0, 'EV': 0, 'EQ': 0, 'EP': 0, 'ES': 0, 'ER': 0, 'II': 0, 'VQ': 0, 'IK': 0, 'VT': 0, 'IN': 0, 'KC': 0, 'KA': 0, 'KG': 0, 'KF': 0, 'KE': 0, 'KD': 0, 'KK': 0, 'KI': 0, 'KH': 0, 'KN': 0, 'KM': 0, 'KL': 0, 'KS': 0, 'KR': 0, 'KQ': 0, 'KP': 0, 'KW': 0, 'KV': 0, 'KT': 0, 'KY': 0, 'DN': 0, 'DL': 0, 'DM': 0, 'DK': 0, 'DH': 0, 'DI': 0, 'DF': 0, 'DG': 0, 'DD': 0, 'DE': 0, 'DC': 0, 'DA': 0, 'DY': 0, 'DV': 0, 'DW': 0, 'DT': 0, 'DR': 0, 'DS': 0, 'DP': 0, 'DQ': 0, 'QQ': 0, 'QP': 0, 'QS': 0, 'QR': 0, 'QT': 0, 'QW': 0, 'QV': 0, 'QY': 0, 'QA': 0, 'QC': 0, 'QE': 0, 'QD': 0, 'QG': 0, 'QF': 0, 'QI': 0, 'QH': 0, 'QK': 0, 'QM': 0, 'QL': 0, 'QN': 0, 'WG': 0, 'WF': 0, 'WE': 0, 'WD': 0, 'WC': 0, 'WA': 0, 'WN': 0, 'WM': 0, 'WL': 0, 'WK': 0, 'WI': 0, 'WH': 0, 'WW': 0, 'WV': 0, 'WT': 0, 'WS': 0, 'WR': 0, 'WQ': 0, 'WP': 0, 'WY': 0, 'PR': 0, 'PS': 0, 'PP': 0, 'PQ': 0, 'PV': 0, 'PW': 0, 'PT': 0, 'PY': 0, 'PC': 0, 'PA': 0, 'PF': 0, 'PG': 0, 'PD': 0, 'PE': 0, 'PK': 0, 'PH': 0, 'PI': 0, 'PN': 0, 'PL': 0, 'PM': 0, 'CK': 0, 'CI': 0, 'CH': 0, 'CN': 0, 'CM': 0, 'CL': 0, 'CC': 0, 'CA': 0, 'CG': 0, 'CF': 0, 'CE': 0, 'CD': 0, 'CY': 0, 'CS': 0, 'CR': 0, 'CQ': 0, 'CP': 0, 'CW': 0, 'CV': 0, 'CT': 0, 'IY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0, 'VF': 0, 'VG': 0, 'IQ': 0, 'VI': 0, 'IS': 0, 'IR': 0, 'VL': 0, 'VM': 0, 'IW': 0, 'IV': 0, 'VP': 0, 'IH': 0, 'VR': 0, 'VS': 0, 'IM': 0, 'IL': 0, 'VV': 0, 'VW': 0, 'IA': 0, 'VY': 0, 'IC': 0, 'IE': 0, 'ID': 0, 'IG': 0, 'IF': 0, 'HY': 0, 'HR': 0, 'HS': 0, 'HP': 0, 'HQ': 0, 'HV': 0, 'HW': 0, 'HT': 0, 'HK': 0, 'HH': 0, 'HI': 0, 'HN': 0, 'HL': 0, 'HM': 0, 'HC': 0, 'HA': 0, 'HF': 0, 'HG': 0, 'HD': 0, 'HE': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0, 'NN': 0, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0, 'NG': 0, 'NY': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'TY': 0, 'TV': 0, 'TW': 0, 'TT': 0, 'TR': 0, 'TS': 0, 'TP': 0, 'TQ': 0, 'TN': 0, 'TL': 0, 'TM': 0, 'TK': 0, 'TH': 0, 'TI': 0, 'TF': 0, 'TG': 0, 'TD': 0, 'TE': 0, 'TC': 0, 'TA': 0, 'AA': 0, 'AC': 0, 'AE': 0, 'AD': 0, 'AG': 0, 'AF': 0, 'AI': 0, 'AH': 0, 'AK': 0, 'AM': 0, 'AL': 0, 'AN': 0, 'AQ': 0, 'AP': 0, 'AS': 0, 'AR': 0, 'AT': 0, 'AW': 0, 'AV': 0, 'AY': 0, 'VK': 0}
	dipeptide_frequency = []

	for i in range(0,len(sequence)-1):
		dipeptides_dict[sequence[i:i+2]] +=1.0

	strain_df = []
	i = 1
	for x,y in sorted(dipeptides_dict.items()):
		strain_df.append((i,y))
		i += 1.0

	dipeptide_frequencies.append(strain_df)	

print sorted(amino_acids)
print sorted(dipeptides_dict)

inputData = dipeptide_frequencies
k=15


centroids = random.sample(inputData, k) #creates the initial centroids from two randomly selected members of the inputData
centroid_members = [] #this list will contain empty lists for each centroid that can be appended with the expression values that fall within that cluster
new_centroids = [] #this list will store the new centroid values during the update step 
while len(centroid_members) < k: #this block adds k number of empy lists
	centroid_members.append([]) 
	new_centroids.append([])
iterations = 0 #this keeps track of the iterations, although it's key role is telling the program after one iteration to 


while (new_centroids != centroids): #this loop can end once the centroids do not change between two iterations 
	if iterations>0: #if not the first iteration, the new_centroids become the centroids variable. This allows the centroids created at the bottom to be compared to the centroids (which will break the loop if they are equal), and set as the centroids to be used if they are not yet equal 
		centroids = new_centroids 
		centroid_members = [] #this rests the centroid numbers list, making sure that the program doesn't retain all past centroid lists. List is created same as above 
		while len(centroid_members) < k:
			centroid_members.append([]) 

	for dataPoint in inputData: #goes through each variance in the ordered variances
		distances = {} #creates a dictionary of distances
		for centroid in centroids: #for each of the centroids calculates the distance and adds it to the distances dictionary with distance as key and the associated centroid as the value
			distance=0
			for i in range(0,len(centroid)):
				distance += (((dataPoint[i][0]-centroid[i][0]) + (dataPoint[i][1]-centroid[i][1]))**2) #adds (x2-x1 + y2-y1)**2 to the distance for all points
			distance = math.sqrt(distance) #takes the squareroot of distance, completing the distance formula
			distances[distance] = centroid #stores the distance for the dataPoint to the centroid in a dictionary with the distance as key and centroid positions as values
		minimum = min(distances) #find the minimum distance in the distances dictionary
		nearest = distances[minimum] #gets the centroid associated with the minimum distance in the dictionary 
		index = centroids.index(nearest) #gets the index for the nearest centroid by looking it up in the list of centroids 
		centroid_members[index].append(dataPoint) #adds the dataPoint to the appropriate centroid list using the index found above
		

		
	new_centroids = [] #this will store new centroids 
	while len(new_centroids) < k: #puts k number of empty lists into the list of new centroids to be appended 
		new_centroids.append([])
	for member in centroid_members: #takes each list from the centroid_members, which each contain values associated with that centroid 
		index = centroid_members.index(member) #this gets the index value, which tells us which number centroid we're working with 
		new = []  #this will store the points for the new centroid 
		x_sum = [] #stores the sums for the x values at CON, H9, and H18
		y_sum = [] #stores the sums for the y values 
		
		for group in member: #for each trio (three expression values) in the centroid list (this represents one gene)
			while len(x_sum) < len(group):
				x_sum.append(0)
				y_sum.append(0)
			for i in range(0,len(group)): #this adds the x and y values from each grouping in the list to the x_sum and y_sum lists 
				x_sum[i] += float(group[i][0]) #this is made into a float so that calculations don't just round to the values 0,1,2 later 
				y_sum[i] += float(group[i][1])
		
		x_avg = []	
		y_avg = []
		for i in range(0,len(x_sum)):
			x_avg.append(x_sum[i]/len(member))
			y_avg.append(y_sum[i]/len(member))
		new = []
		for i in range(0,len(x_avg)):
			new.append((x_avg[i],y_avg[i]))
		new_centroids[index] = new #this uses the index obtained earlier to add the new centroid values to the new_centroids list at the correct position 
	iterations +=1 #adds 1 to the current iteration count
	
centroids = new_centroids  #now that the loop is finished, the new_centroids are our final set of centroids 
centroid_members = [] #this establishes the centroid_members one last time 
while len(centroid_members) < k:
	centroid_members.append([])	

for dataPoint in inputData: #Finds the final centroid associations after the while loop obtained the set centroid positions. This code is commented above
	distances = {}
	for centroid in centroids:
		distance=0
		for i in range(0,len(centroid)):
			distance += (((dataPoint[i][0]-centroid[i][0]) + (dataPoint[i][1]-centroid[i][1]))**2) #adds (x2-x1 + y2-y1)**2 to the distance for all points
		distance = math.sqrt(distance)
		distances[distance] = centroid
	minimum = min(distances)
	nearest = distances[minimum]
	index = centroids.index(nearest)
	centroid_members[index].append(dataPoint)
	

for centroid in centroid_members: #gets each centroid group
	for member in centroid: #gets each gene's expression values for the centroid group
		x_values = []
		y_values = []
		for i in range(0,len(member)): #this stores the x and y values for each gene's expression values 
			x_values.append(member[i][0])
			y_values.append(member[i][1])
		plots.plot(x_values,y_values) #this plots all the lines for each centroid 
	plots.show() #this shows the plots

