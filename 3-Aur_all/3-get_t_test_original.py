
import itertools
import os

from scipy.stats import ttest_ind, ttest_ind_from_stats
import pandas as pd
from scipy.stats.stats import pearsonr
import scipy.stats as stats
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
import math
def get_snp_fre():

    f1 = open('/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/Min_data/chang_paper/Aureus_final_new_version/NC_007795.1/result_new_version_cutoff_score3_3_4.5/n1_nodes_vector_original.bed','r')
    lines1 = f1.readlines()

    dict1 = {}

    for line1 in lines1:

        name1 = line1.split('\t')[0].split('_')[0] + ',' + line1.split('\t')[0].split('_')[1] # A,location

        fre_list = [int(i) for i in line1.split('\t')[1].strip()[1:-1].split(', ')]

        dict1[name1] = fre_list

    return dict1


def getsta(list1):

    b = pd.Series(list1)
    c = b.describe()
    print(list(c))

    return list(c)

def get_strain_SNPs(file1,output_path):

    f1 = open(file1,'r')

    lines1 = f1.readlines()

    SNP_fre = get_snp_fre()

    dict1 = {}

    name1 = lines1[0][1:]
    fre_list = []
    for line1 in lines1[1:]:

        fre_list.append(line1.strip())

        # fre_list = sorted([i.split('_')[0] for i in line1.split('\n')[:-1]], key = lambda x:int(x.split('_')[1]), reverse= True)
    fre_list = sorted([i for i in fre_list], key=lambda x: int(x.split('_')[1]),reverse=True)
    print(fre_list)
    fre_list1 = []

    for i in fre_list:
        if i.split('_')[0] in list(SNP_fre.keys()):
            fre_list1.append(i)

    dict1[name1] = fre_list1

    columns1 = ['SNP1', 'SNP2', 'r', 'P-value']




    for key,val in dict1.items():
        f1 = open(output_path + '/' + key,'w')
        f1.write('\t'.join(columns1) + '\n')
        for name1, name2 in itertools.combinations(val[:1000], 2):
            print(name1,name2)
            fre1 = SNP_fre[name1.split('_')[0]]
            fre2 = SNP_fre[name2.split('_')[0]]

            r = pearsonr(fre1,fre2)[0]

            if r !=1:
                temp_value = r*math.sqrt(66/(1-r*r))

                p = 1 - stats.t.cdf(temp_value, 66)

                f1.write(name1 + '\t' + name2 + '\t' + str(round(r,6)) + '\t' + str(round(p,6)) + '\n')
            else:
                f1.write(name1 + '\t' + name2 + '\t' + str(round(r,6)) + '\t' + 'NA' + '\n')

        f1.close()



input_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/known_strain_unique_SNP/Aur'
output_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/t-test_original/Aur'
files = [i for i in os.listdir(input_path) if 'SNPs_statistics.bed' not in i]

for file in files:

        input_file = input_path + '/' + file

        output_path1 = output_path

        if os.path.exists(output_path1) == False:
            os.makedirs(output_path1)

        get_strain_SNPs(input_file,output_path1)










