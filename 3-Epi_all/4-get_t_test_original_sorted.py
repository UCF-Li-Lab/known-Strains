
import itertools
import os

from scipy.stats import ttest_ind, ttest_ind_from_stats
import pandas as pd
from scipy.stats.stats import pearsonr
import scipy.stats as stats
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
import math

def get_sorted(file1,file2):

    f1 = open(file1,'r')
    f2 = open(file2,'w')

    lines1 = f1.readlines()
    all_list = []
    columns1 = ['SNP1', 'SNP2', 'r', 'P-value']
    f2.write('\t'.join(columns1) + '\n')

    for line1 in lines1[1:]:

        all_list.append(line1.strip().split('\t'))


    all_list = sorted(all_list, key =lambda x:float(x[2]),reverse=True)

    for i in all_list:
        f2.write('\t'.join(i) + '\n')

    f1.close()
    f2.close()











input_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/t-test_original/Epi'
output_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/t-test_original_sorted/Epi'


files = os.listdir(input_path)

for file in files:
    input_file = input_path + '/' + file
    output_file = output_path + '/' + file

    get_sorted(input_file,output_file)










