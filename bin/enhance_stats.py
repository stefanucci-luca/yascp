#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np

from scipy import stats

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-id", "--donor_id", metavar="donor_id", dest = "donor_id",
                    help="", type=str,
                    default=1)
parser.add_argument("-gp", "--genotype_phenotype_mapping", metavar="genotype_phenotype_mapping", dest = "genotype_phenotype_mapping",
                    help="", type=str,
                    default=None)

parser.add_argument("-if", "--input_file", metavar="input_file", dest = "input_file",
                    help="", type=str,
                    default=None)

parser.add_argument("-dm", "--dm", metavar="input_file", dest = "donor_match",
                    help="", type=str,
                    default=None)

parser.add_argument("-m", "--mode", metavar="input_file", dest = "mode",
                    help="", type=str,
                    default=None)
                    
args = parser.parse_args()
mode =args.mode
if (args.genotype_phenotype_mapping):
    gp_ma = pd.read_csv(args.genotype_phenotype_mapping,sep='\t',index_col=0)
    gp_ma.index = gp_ma.index.astype(str)
else:
    gp_ma = pd.DataFrame()
input_id =args.donor_id
# input_id = 'ELGH_VAL11650907'
# stats_CRD_CMB12979963_gt_donor_assignments.csv
GT_Assignments = pd.read_csv(args.donor_match,sep=',')
GT_Assignments=GT_Assignments.set_index('donor_query')
import os
tranche = os.getcwd().split('/')[-4]
pool = input_id
input_table_file = pd.read_csv(args.input_file,sep='\t')
input_table_file =input_table_file.set_index('experiment_id')


Unassigned= []
try:
    All_expected_ids = input_table_file.loc[args.donor_id,'donor_vcf_ids'].replace('\'','').split(',')
except:
    All_expected_ids =[]
Good_ids=[]
Emergency_ids=[]

for id1 in All_expected_ids:
    if '-999-' in id1:
        Emergency_ids.append(id1)
    else:
        Good_ids.append(id1)
Good_ids = ",".join(Good_ids)
Emergency_ids = ",".join(Emergency_ids)
GT_Assignments['Emergency_ids expected']=Emergency_ids
GT_Assignments['Good_ids expected']=Good_ids
GT_Assignments['Match Expected']='False'
GT_Assignments['pool']=pool
GT_Assignments['tranche']=tranche
GT_Assignments['donor_gt original']=GT_Assignments['donor_gt']
D2 = pd.DataFrame(All_expected_ids,columns=['col1'])
All_Expected_set = []
for ix in GT_Assignments.index:
    # print(ix)
    #Here should add a a filter to estimate whether it is a good match.
    replacement = GT_Assignments.loc[ix,'donor_gt']
    # if (replacement=='15001506223417_203765330119_R02C02'):
    #     print('check')
    expected = 'NA'
    poor_replacement =''
    if ix in Unassigned:
       poor_replacement = 'Poor_GT_score___'

    # The folowing will attempt to replcat the donor names as per phenotype file, as provided by mapping.
    # If it fails it will keep the genotype id.
    
    if (len(gp_ma)>0):
        try:
            # replacement = 'S2-046-01741'
            try:
                replacement = str(gp_ma.loc[replacement].iloc[0].values[0])
            except:
                replacement = str(gp_ma.loc[replacement].iloc[0])
            # print(replacement)
            # replacement = replacement.values[0]
            
            if len(D2[D2.col1.str.contains(replacement)])>0:
                GT_Assignments.loc[ix,'Match Expected']='True'
                remove_from_set = D2[D2.col1.str.contains(replacement)]['col1'].values[0]
                All_Expected_set.append(remove_from_set)
        except:
            try:
                # replacement='15001608190388_204238910153_R09C02'
                
                replacements = gp_ma.loc[replacement.split('_')[0]]
                
                if len(replacements)>1:
                    for rep1 in replacements.iloc[:,0]:
                        
                        if len(D2[D2.col1.str.contains(rep1)])>0:
                            GT_Assignments.loc[ix,'Match Expected']='True'
                            replacement = rep1
                            remove_from_set = D2[D2.col1.str.contains(replacement)]['col1'].values[0]
                            All_Expected_set.append(remove_from_set)
                else:
                    replacement = replacements.values[0]
            except:
                replacement = 'No_mapping___'+replacement
    else:
        _ = ''
    if len(D2[D2.col1.str.contains(replacement)])>0:
        GT_Assignments.loc[ix,'Match Expected']='True'
        remove_from_set = D2[D2.col1.str.contains(replacement)]['col1'].values[0]
        All_Expected_set.append(remove_from_set)
    if 'THP1' in replacement:
        replacement = 'celline_THP1'
    if 'U937' in replacement:
        replacement = 'celline_U937'
    if poor_replacement == 'Poor_GT_score___':
        replacement = poor_replacement+replacement

    GT_Assignments.loc[ix,'donor_gt']=replacement


# Now that we have loaded them replace all the files with the correct donor ids.
Missing = ';'.join(set(All_expected_ids)-set(All_Expected_set))
GT_Assignments['Missing']=Missing
GT_Assignments.to_csv(f'GT_replace_{args.donor_match}',sep='\t')



print("Done")