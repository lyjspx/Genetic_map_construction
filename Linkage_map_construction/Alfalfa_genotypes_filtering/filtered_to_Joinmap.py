input_file = 'filter_parent_polymor_output.txt'
outpuit_file = 'data_for_Joinmap.txt'
write_buffer = ''
line_list = list()
with open(input_file) as input_filef:
	write_buffer += input_filef.readline()
	for line in input_filef:
		list_temp = list()
		aline = line.strip().split('\t')
		list_temp.append(aline[0])
		if aline[1] not in ['A','G','C','T']  and aline[2] not in ['A','G','C','T']:
			print(aline[0])
		elif aline[1] not in ['A','G','C','T'] and aline[2] in ['A','G','C','T']:
			for snp in aline[3:]:
				if snp == 'N':
					list_temp.append('-')
				elif snp not in ['A','G','C','T']:
					list_temp.append('h')
				elif snp == aline[2]:
					list_temp.append('b')
				elif snp in ['A','G','C','T']:
					list_temp.append('a')
				else:
					print(aline[0])
			list_temp.insert(1,'b') #insert parent_b genotype
			list_temp.insert(1,'h')
		elif aline[1] in ['A','G','C','T'] and aline[2] not in ['A','G','C','T']:
			for snp in aline[3:]:
				if snp == 'N':
					list_temp.append('-')
				elif snp not in ['A','G','C','T']:
					list_temp.append('h')
				elif snp == aline[1]:
					list_temp.append('a')
				elif snp in ['A','G','C','T']:
					list_temp.append('b')
				else:
					print(aline[0])
			list_temp.insert(1,'h') #insert parent_b genotype
			list_temp.insert(1,'a')
		elif aline[1] in ['A','G','C','T']and aline[2] in ['A','G','C','T']:
			for snp in aline[3:]:
				if snp == 'N':
					list_temp.append('-')
				elif snp not in ['A','G','C','T']:
					list_temp.append('h')
				elif snp == aline[1]:
					list_temp.append('a')
				elif snp == aline[2]:
					list_temp.append('b')
				else:
					print(aline[0])
			list_temp.insert(1,'b') #insert parent_b genotype
			list_temp.insert(1,'a')

		else:
			print(aline[0])
		line_list.append(list_temp)
		#print(aline)
for oneline in line_list:
	write_buffer += oneline[0]
	write_buffer += '\t'
	for each in oneline[1:]:
		write_buffer += each.lower()
		write_buffer += '\t'
	write_buffer += '\n'
with open(outpuit_file, 'w') as outpuit_filef:
	outpuit_filef.write(write_buffer)
print('Finished')

		
		






