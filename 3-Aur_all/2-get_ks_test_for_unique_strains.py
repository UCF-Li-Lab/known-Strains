import os

import os

def get_strain_SNPs(file1):

    f1 = open(file1,'r')

    lines1 = f1.readlines()

    lines1 = ''.join(lines1).split('>')[1:]

    dict1 = {}
    for line1 in lines1:

        name1 = line1.split('\n')[0]
        # print(line1.split('\n')[:-1])
        fre_list = [int(i.split('_')[1]) for i in line1.split('\n')[1:-1]]

        dict1[name1] = fre_list

    return dict1




def get_ks(list1):
    import scipy.stats as ss
    A = list1
    a = ss.kstest(A,'uniform', args=(min(A),max(A)-min(A)))
    b = ss.kstest(list(set(A)),'uniform', args=(min(A),max(A)-min(A)))
    print(list(a))
    print(list(b))

    return [list(a)[1],list(b)[1]]



input_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/known_strain_unique_SNP/Aur'
output_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/ks_test/Aur'

files = [i for i in os.listdir(input_path) if 'SNPs_statistics.bed' not in i]


f3 = open(output_path + '/' + 'Aur_original.txt','w')
f4 = open(output_path + '/' + 'Aur_non_repeat.txt','w')
for file in files:
    input_file = input_path + '/' + file
    output_file = output_path + '/' + file.split('.txt')[0] + '_original.txt'
    output_file1 = output_path + '/' + file.split('.txt')[0] + '_non_repeat.txt'
    dict1 = get_strain_SNPs(input_file)
    print(file)

    f1 = open(output_file,'w')
    f2 = open(output_file1, 'w')
    for key,val in dict1.items():
        print(key)
        list1 = get_ks(val)
        f1.write(key + '\t' + str(list1[0]) + '\n')
        f2.write(key + '\t' + str(list1[1]) + '\n')

        f3.write(key + '\t' + str(list1[0]) + '\n')
        f4.write(key + '\t' + str(list1[1]) + '\n')

    f1.close()
    f2.close()


f3.close()
f4.close()