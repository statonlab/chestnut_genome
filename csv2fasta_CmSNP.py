




# fileName = 'SNP_MOESM2_ESM.csv'
# fileName = 'SSR_MOESM2_ESM.csv'

def snp(seq,start = 0):
	ind_0 = seq.index('[')
	part1 = seq[0:ind_0]
	part2 = seq[ind_0+5:]
	if start == 0:
		print "I am selecting the first base......"
		middle_part = seq[ind_0+1]
	else:
		print "I am selecting the second base......"
		middle_part = seq[ind_0+3]
	new_seq = part1 + middle_part + part2
	return new_seq

def generateFASTA(fasta_content):
	new_file_name = fileName + '.1.' +'fasta'
	with open(new_file_name, 'w') as fastaFile:
		for i in fasta_content:
			new_string  = '\n'
			new_string = new_string + '>' + i[0]
			new_string = new_string + '\n'
			new_string = new_string + i[1] + '\n'
			fastaFile.write(new_string)
		# print new_string



def main():
	with open(fileName) as file:
		content = file.readline()
	content_matrix = content.split('\r')
	fasta_content = []
	header = 0
	for i in content_matrix:
		if header == 0:
			header = header + 1
			continue
		# print "\n"
		name = i.split(',')[0]
		seq = i.split(',')[-2]
		new_seq = snp(seq,start=0)
		fasta_content.append([name,new_seq])
		# print name, new_seq
	generateFASTA(fasta_content)












if __name__ == '__main__':
	main()




















