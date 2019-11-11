vcf_file_name = 'rawSNP.vcf'
output_hapmap_count_file = 'Hapmap.hmc.txt'
count_container = list()
SNP_name_container = list()
sample_name = list()
final_list = list()
with open(vcf_file_name) as vcf:
        for line in vcf:
                line_count = list()
                if line[0:2] == '#C': sample_name = line.strip().split('\t')[9:]
                elif line[0] != '#':
                        line_list = line.strip().split('\t')
                        SNP_name_container.append(line_list[2])
                        line_list = line_list[9:]
                        #print(line_list)
                        for element in line_list:
                                #print(element)
                                if len(element.split(':')) == 1:
                                    allele1, allele2 = 0, 0
                                elif len(element.split(':')[1].split(',')) == 2:       
                                        allele1, allele2 = element.split(':')[1].split(',')[0], element.split(':')[1].split(',')[1]
                                else:
                                        allele1, allele2 = element.split(':')[1].split(',')[0], 0
                                temp = [allele1,'|',allele2]
                                line_count.append(temp)
                        count_container.append(line_count)
#sample_name[-1]=sample_name[:-2]
sample_name.append('HettCount_allele1')
sample_name.append('HettCount_allele2')
sample_name.append('Count_allele1')
sample_name.append('Count_allele2')
sample_name.append('Frequency')
for line_count_info in count_container:
        final_line = list()
        heter_count_allele1, heter_count_allele2, count_allele1, count_allele2 = 0, 0, 0, 0        
        for oneSNP in line_count_info:
                referSNP, altSNP = int(oneSNP[0]), int(oneSNP[2])
                if referSNP > 0 and altSNP > 0: 
                        heter_count_allele1 += referSNP
                        heter_count_allele2 += altSNP
                count_allele1 += referSNP
                count_allele2 += altSNP
        final_line += line_count_info
        final_line.append(heter_count_allele1)
        final_line.append(heter_count_allele2)
        final_line.append(count_allele1)
        final_line.append(count_allele2)
        if count_allele1 + count_allele2 != 0:
            final_line.append(count_allele1/(count_allele2+count_allele1))
        else: final_line.append(0)
        final_list.append(final_line)
        
with open(output_hapmap_count_file, 'w') as hmc:
        hmc.write('rs')
        for x in sample_name:
                hmc.write('\t' + x)
        hmc.write('\n')
        for SNP in SNP_name_container:
                hmc.write(SNP)
                for count in final_list[SNP_name_container.index(SNP)]:
                    if type(count) is list:    
                        hmc.write('\t' + str(count[0])+str(count[1])+str(count[2]))
                    else:hmc.write('\t' + str(count))
                hmc.write('\n')

                

