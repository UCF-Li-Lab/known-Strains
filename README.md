# KnownStrainPaper_CodesUsed

The following details are an explanation of the in-house codes used for analysis in the paper.
These codes were used in the process of analyzing the presence of the known strains previoudly identified
by Pathoscope2 and StrainEst in the metagenomics dataset.
 
Step-1:
This code was used to get the SNPs of the known strains identified Pathoscope 2 and StrainEst tools.

4-known_strain_mutation

Requirements:
(a) genome sequence of the species of interest (example: S. aureus or S. epidermidis)
(b) genome sequences of the known strains identified

Output: SNPs of the known strains 

Step 2:
comprises the analysis of unique SNPs of known strains, their coverages and statistics.
The following set of codes are included in each of the folder labeled as: 3-Aur_all and 3-Epi_all

Step-2_A:
1-get_unique_SNPs
This was used to get the unique SNPs of known strains.
Unique SNPs here means only occurs in one known strains and not shared.

Requirement:
(a) the ouput from Step-1 (known strain SNPs)

Output: Unique SNPs of known strain and their coverages

Step-2_B:
2-get_ks_test_for_unique_strains
This code was used to determine uniform distribution using KS test

Requirement: output from Step-2_A (Unique SNPs of known strain and their coverages)
Output: KS-test results

Step-2_C:
3-get_t_test_original
This code was used to perform t-test analysis

Requirement: output from Step-2_A (Unique SNPs of known strain and their coverages)
Output: t-test results

Step-2_D:
4-get_t_test_original_sorted
This code was used to sort the t-test results from Step-2_C

Requirement: output from Step-2_C (t-test results)
Output: sorted t-test results

Step-2_E:
5-get_get_significant_p_value_nums
This code was used to test significance

Requirement: output from Step-2_D (sorted t-test)
Output: significance result
