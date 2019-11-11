input_file_name = 'filter_missing_parent_output.txt'
out_file_name = 'filter_parent_polymor_output.txt'
line_container = list()
with open(input_file_name) as input_file_name_f:
	line_container.append(input_file_name_f.readline())
	for line in input_file_name_f:
		one_line = list()
		one_line = line.split('\t')
		if one_line[1] != one_line[2]: 
			line_container.append(line)

with open(out_file_name, 'w') as out_file_name_f:
	for line in line_container:
		out_file_name_f.write(line) 