input_hap_count_file = 'Hapmap.hmc.txt'
input_filtered_file = 'filter_parent_polymor_output.txt'
output_hap_count_file = 'filterd_hap_count.txt'

name_list = list()
write_buffer = ''
with open(input_filtered_file) as filered:
	next(filered)
	for line in filered:
		name_list.append(line.strip().split('\t')[0])

with open(input_hap_count_file) as hap_count:
	write_buffer += hap_count.readline()
	for line in hap_count:
		if line.strip().split('\t')[0] in name_list:
			write_buffer += line

with open(output_hap_count_file, 'w') as output:
	output.write(write_buffer)
