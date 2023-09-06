
import os
import pandas as pd
def get_number(file1):
    f1 = open(file1,'r')
    lines1 = f1.readlines()

    num_all =len(lines1)-1

    num1 = 0
    p_all = 0
    for line1 in lines1[1:]:
        p_value = line1.split('\t')[-1].strip()

        if p_value == 'NA':
            p_value1 = 0
            if p_all < 0.01:
                p_all += p_value1
                num1 +=1
            else:
                if p_all >= 0.01:
                    break
        else:
            p_value1 = float(p_value)
            if p_all < 0.01:
                p_all += p_value1
                num1 +=1
            else:
                if p_all >= 0.01:
                    break

    num2 = num1-1

    per = round(num2/num_all,4)

    print(num_all,num2,per)

    return [num_all,num2,per]






input_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/t-test_original_sorted/Aur'

files = os.listdir(input_path)

data = []
index1 = []
collumns1 = []

for file in files:
    print(file)
    index1.append(file)
    input_file = input_path + '/' + file
    list1 = get_number(input_file)

    data.append(list1)

collumns1 = ['#pairs','#significant pairs','percentage of significant pairs']

df = pd.DataFrame(data,index=index1,columns=collumns1)

df.to_csv('/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/t-test_original_sorted/Aur_stat',sep = '\t')






