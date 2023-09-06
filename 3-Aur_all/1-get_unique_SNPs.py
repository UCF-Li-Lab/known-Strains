


import os
import pandas as pd

def read_file(file1):
    f1 = open(file1,'r')
    lines1 = f1.readlines()

    dict1 = {}

    for line1 in lines1:
        loc = line1.split('\t')[0]
        A_fre = line1.split('\t')[1]
        C_fre = line1.split('\t')[2]
        G_fre = line1.split('\t')[3]
        T_fre = line1.split('\t')[4].strip()

        if int(A_fre) != 0:
            dict1['A' + ',' + loc] = int(A_fre)
        if int(C_fre) != 0:
            dict1['C' + ',' + loc] = int(C_fre)
        if int(G_fre) != 0:
            dict1['G' + ',' + loc] = int(G_fre)
        if int(T_fre) != 0:
            dict1['T' + ',' + loc] = int(T_fre)

    return dict1







def get_predict_SNPdict(input_path,filer_file,output_path):
    files = [i for i in os.listdir(input_path) if 'txt' not in i]

    dict_orginal_SNP = {}
    for file in files:
        input_file = input_path + '/' + file + '/mutation_file'

        f1 = open(input_file,'r')
        lines1 = f1.readlines()
        lines1 = ''.join(lines1).split('>')[1:]
        for line1 in lines1:
            strain = line1.split('\n')[0]
            SNPs = [i for i in line1.split('\n')[1:-1] if '.' not in i]
            dict_orginal_SNP[file + '--' + strain] = list(set(SNPs))

    all_SNPs = []

    for key, val in dict_orginal_SNP.items():
        if key != 'Pathoscope2AurKS_genomes--CIGC93_GCF_000249035.1_ASM24903v2_genomic.fasta' and key !='Pathoscope2AurKS_genomes--USA300TCH959_GCF_000153665.1_ASM15366v1_genomic.fasta' \
                and key != '' and key!= '':
            all_SNPs += val

    from collections import Counter

    result = Counter(all_SNPs)

    unique_SNPs = []
    # two_SNPs = []
    for key, val in result.items():

        if val == 1:
            unique_SNPs.append(key)

        # if val == 2:
        #     two_SNPs.append(key)
    input_file1 = filer_file
    fre_dict = read_file(input_file1)

    dict1 = {}
    data = []
    index1 = []

    for key,val in dict_orginal_SNP.items():
        strain = key
        SNPs = [i for i in val if i in unique_SNPs]
        if len(SNPs) != 0 :
            SNPs_fre = [fre_dict[i] for i in SNPs if i in fre_dict.keys() and fre_dict[i]!=0]
            aver = round(sum(SNPs_fre)/len(SNPs_fre),4)

            selected_SNPs = [i + '_' + str(fre_dict[i]) for i in SNPs if i in fre_dict.keys() and fre_dict[i]!=0] # format A,location_frequency
            dict1[strain] = [str(aver),selected_SNPs]

            index1.append(strain)
            data.append([len(val),len(SNPs),len(SNPs_fre),len(SNPs)-len(SNPs_fre)])

            f2 = open(output_path + '/' + strain,'w')

            f2.write('>' + strain + '--' + str(aver) + '\n')  ## it's the original average
            f2.write('\n'.join(selected_SNPs) + '\n')

        else:
            # SNPs_fre = [fre_dict[i] for i in SNPs if i in fre_dict.keys() and fre_dict[i] != 0]
            # aver = round(sum(SNPs_fre) / len(SNPs_fre), 4)

            # selected_SNPs = [i + '_' + str(fre_dict[i]) for i in SNPs if
            #                  i in fre_dict.keys() and fre_dict[i] != 0]  # format A,location_frequency
            # dict1[strain] = [str(aver), selected_SNPs]

            index1.append(strain)
            data.append([len(val), len(SNPs), 0, 0])

            f2 = open(output_path + '/' + strain, 'w')

            f2.write('>' + strain + '--' + 'NA' + '\n')  ## it's the original average
            # f2.write('\n'.join(selected_SNPs) + '\n')

    stat_file = output_path + '/SNPs_statistics.bed'

    df = pd.DataFrame(data,index=index1,columns=['#Original SNP','#unique SNP','#unique SNP within reads','#unique SNP without reads'])

    df.to_csv(stat_file,sep='\t')

    # dict1 = {}
    # for file in files:
    #     input_file = input_path + '/' + file + '/mutation_file'
    #     input_file1 = filer_file
    #     output_file = output_path + '/' + file + '_unique_SNP.txt'
    #     fre_dict = read_file(input_file1)
    #
    #     f1 = open(input_file,'r')
    #     lines1 = f1.readlines()
    #     lines1 = ''.join(lines1).split('>')[1:]
    #
    #
    #     temp_dict = {}
    #
    #     for line1 in lines1:
    #         strain = line1.split('\n')[0]
    #         SNPs = line1.split('\n')[1:-1]
    #         temp_dict[strain] = SNPs
    #
    #     all_SNPs = []
    #
    #     for key,val in temp_dict.items():
    #         all_SNPs += val
    #
    #
    #
    #     from collections import Counter
    #
    #     result = Counter(all_SNPs)
    #
    #     unique_SNPs = []
    #     two_SNPs = []
    #     for key,val in result.items():
    #
    #         if val ==1:
    #             unique_SNPs.append(key)
    #
    #         if val ==2:
    #             two_SNPs.append(key)
    #
    #     f2= open(output_file,'w')
    #
    #     for key,val in temp_dict.items():
    #         strain = key
    #         SNPs_test = [i for i in val if i in unique_SNPs]
    #         if len(SNPs_test) != 0:
    #             SNPs = [i for i in val if i in unique_SNPs]
    #             SNPs_fre = [fre_dict[i] for i in SNPs if i in fre_dict.keys() and fre_dict[i]!=0]
    #             aver = round(sum(SNPs_fre)/len(SNPs_fre),4)
    #
    #             selected_SNPs = [i + '_' + str(fre_dict[i]) for i in SNPs if i in fre_dict.keys() and fre_dict[i]!=0] # format A,location_frequency
    #             dict1[strain] = [str(aver),selected_SNPs]
    #
    #             f2.write('>' + strain + '--' + str(aver) + '\n')
    #             f2.write('\n'.join(selected_SNPs) + '\n')
    #
    #         else:
    #             SNPs = [i for i in val if i in two_SNPs]
    #             SNPs_fre = [fre_dict[i] for i in SNPs if i in fre_dict.keys() and fre_dict[i]!=0]
    #             aver = round(sum(SNPs_fre)/len(SNPs_fre),4)
    #
    #             selected_SNPs = [i + '_' + str(fre_dict[i]) for i in SNPs if i in fre_dict.keys() and fre_dict[i]!=0] # format A,location_frequency
    #             dict1[strain] = [str(aver),selected_SNPs]
    #             f2.write('>' + strain + '--' + str(aver) + '\n')
    #             f2.write('\n'.join(selected_SNPs) + '\n')
    #
    #     f1.close()
    #     f2.close()
    # return dict1



input_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/known_strain_SNPs/Aur'
filter_file = '/media/saidi/Elements/Project/Project17_1_MSMS/testing_real/Min_data/chang_paper/Aureus_final_new_version/NC_007795.1/merged_filter_polymorphic_sites_100_filter'
output_path = '/media/saidi/Elements/Project/Project_for_Min/KSPaper_by_me1/known_strain_unique_SNP/Aur'


predict_strain = get_predict_SNPdict(input_path,filter_file,output_path)













