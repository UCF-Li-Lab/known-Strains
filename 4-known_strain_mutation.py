import argparse
import os, sys
# Initialize parser

'''
python known_strain_mutation.py --ref ref_file --known_strain_path known_path --output  output_path
'''


parser = argparse.ArgumentParser()

parser.add_argument('--ref', help='Give a unique sample name',
                    default=None)
parser.add_argument('--known_strain_path', help='Give a unique sample name',
                    default=None)
parser.add_argument('--output', help='all the comparison result',
                    default=True)

args = parser.parse_args()

# ref = '/media/saidi/Elements/Project/Project17_mixtureS_from_Xin_orginal/latest/MixtureS/Ming/known_strain_mutation/test_data/ref/ref.fna'
# known_strain_path = '/media/saidi/Elements/Project/Project17_mixtureS_from_Xin_orginal/latest/MixtureS/Ming/known_strain_mutation/test_data/known_strains'
# output_path = '/media/saidi/Elements/Project/Project17_mixtureS_from_Xin_orginal/latest/MixtureS/Ming/known_strain_mutation/test_result'


ref = args.ref
known_strain_path = args.known_strain_path
output_path = args.output

print('ref_file:')
print(ref)
print('known_strain_path:')
print(known_strain_path)
print('output path')
print(output_path)



############### get the known strain file location#################
def getKnownStrain(known_strain_path):
    loc_list = []
    files = os.listdir(known_strain_path)
    for file in files:
        loc_list.append(known_strain_path + '/' + file)
    return loc_list

############ align strain sequence with reference ################
def Align_seq(loc_list,ref,output_path):
    for loc in loc_list:
        strain_file = loc
        ref_file = ref
        output_path1 = output_path + '/align_result'
        if os.path.exists(output_path1) == False:
            os.makedirs(output_path1)

        output_path2 = output_path1 + '/' + loc.split('/')[-1]
        if os.path.exists(output_path2) == False:
            os.makedirs(output_path2)

        output_file = output_path2 + '/' + loc.split('/')[-1]

        command = 'dnadiff -p %s %s %s' % (output_file,ref_file,strain_file)

        os.system(command)

############# get mutation ##################

def getMutation(SNP_file):
    f1 = open(SNP_file,'r')
    lines1 = f1.readlines()
    SNP_list = []
    for line1 in lines1:
        loc = int(line1.split('\t')[0]) - 1
        snp = line1.split('\t')[2]
        SNP_list.append(','.join([snp,str(loc)]))

    return SNP_list

########## do comparison ################

print('##################### getting the known strain location ###########################')
loc_list = getKnownStrain(known_strain_path)

print('##################### align known strain to reference ###########################')
Align_seq(loc_list,ref,output_path)

SNP_path = output_path + '/align_result'

snp_files = os.listdir(SNP_path)

f1 = open(output_path + '/mutation_file','w')
for file in snp_files:
    snp_file = SNP_path + '/' + file + '/' + file + '.snps'
    snp_list = list(set(getMutation(snp_file)))
    f1.write('>' + file + '\n' + '\n'.join(snp_list) + '\n')

f1.close()



