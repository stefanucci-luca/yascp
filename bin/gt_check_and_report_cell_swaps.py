#!/usr/bin/env python

import pandas as pd
import os
import glob

# TYhis code takes the vireo thats run with 100% of SNPs and itterates through each of the subsampled SNP panels to 
# determine the cells that have swapped the donor identity. Based on the iterations it estimates the error rate for the cell to be 
# assigned to the wrong donor.

all_vireo_runs = glob.glob('./vireo_*/')
# Reference vireo
# All other vireo runs
all_vireo_runs = pd.DataFrame(all_vireo_runs,columns=['col'])
all_vireo_runs['name']=all_vireo_runs['col'].str.split('/').str[-2]
all_vireo_runs['rand itteration']=all_vireo_runs['name'].str.split('___').str[1]
ref_vireo_run =all_vireo_runs[all_vireo_runs['rand itteration'].isna()]
subsampling_itterations = all_vireo_runs[~all_vireo_runs['rand itteration'].isna()]

# Now that we know which are the subsampling itterations and which is the 100% snp itteration we check the samples that swap the identities.
# Since users may be running pipeline in genotype absent mode, we should perform a GT match for every pair to determine which is the matching donor for each.
reference_vcf = f"{list(ref_vireo_run['col'])[0]}GT_donors.vireo.vcf.gz"
out_dir = '/'.join(reference_vcf.split('/')[:-1])
os.system(f"bash fix_vireo_header.sh {reference_vcf} {out_dir}")
reference_vcf = f"{list(ref_vireo_run['col'])[0]}headfix_srt_vcf.vcf.gz"

#Load the reference cell identity file.
Reference_cell_identities = pd.read_csv(f"{list(ref_vireo_run['col'])[0]}donor_ids.tsv",sep='\t')
Reference_cell_identities2=Reference_cell_identities.set_index('cell')
Reference_cell_identities2['Nr times becoming Doublet in subsampling']=0
Reference_cell_identities2['Nr times becoming Unassigned in subsampling']=0
Reference_cell_identities2['Nr times becoming different donor in subsampling']=0
Reference_cell_identities2['New Donor Identities'] = ''
Reference_cell_identities2['Nr vireo subsampling itterations']=0

for subsampling_file in subsampling_itterations['col']:
    print(subsampling_file)
    subsampling_vcf=f"{subsampling_file}GT_donors.vireo.vcf.gz"
    Cell_identities_subsampled = pd.read_csv(f"{subsampling_file}donor_ids.tsv",sep='\t')

    # Perform gtcheck
    out_dir = '/'.join(subsampling_vcf.split('/')[:-1])
    os.system(f"bash fix_vireo_header.sh {subsampling_vcf} {out_dir}")
    subsampling_vcf = f"{subsampling_file}headfix_srt_vcf.vcf.gz"
    os.system(f'bcftools gtcheck --no-HWE-prob -g {reference_vcf} {subsampling_vcf} > {subsampling_file}gt_score_check.tsv')
    # Now that we have performed the GT check we look at how much the cell identities have changed for each of the donors between the refeence and the 
    # subsampled datasets for each of the cells. 
    GT_Data = pd.read_csv(f'{subsampling_file}gt_score_check.tsv',sep='\t',comment='#', skiprows=15)
    # Now loop through each of the donors and see which have changed identities.
    GT_Data.columns=['DC','[2]Query Sample','[3]Genotyped Sample','[4]Discordance','[5]-log P(HWE)','[6]Number of sites compared']
    Cell_identities_subsampled2=Cell_identities_subsampled.set_index('cell')
    for donor1 in set(GT_Data['[2]Query Sample']):
        don_data = GT_Data[GT_Data['[2]Query Sample']==donor1]
        # Need to drop these cells from output if the donor id contains donor.
        if 'donor' in donor1:
            # If we do perform a posthoc gt check then this would need be ok.
            _='this donor was not in bridging file and hence was not exected'
            _='for this reason instead of performing posthoc gt check we drop these from subsequent analysis as these wouldnt be returned anyways.'
            Reference_cell_identities2 = Reference_cell_identities2[Reference_cell_identities2['donor_id']!=donor1]
            continue
        Best_match = list(don_data[don_data['[4]Discordance'] == min(don_data['[4]Discordance'])]['[2]Query Sample'])[0]
        
        
        #Now that we have a best match and the donor of interest we can see how many cells have changed the donor identities.
        Donor_Original_Cells = set(Reference_cell_identities[Reference_cell_identities['donor_id']==donor1]['cell'])
        Donor_Subsampled_Cells = set(Cell_identities_subsampled[Cell_identities_subsampled['donor_id']==Best_match]['cell'])
        Donor_Subsampled_Cell_df = Cell_identities_subsampled[Cell_identities_subsampled['donor_id']==Best_match]
        
        Donor_Subsampled_Cell_df = Donor_Subsampled_Cell_df.set_index('cell')
        Missing_Cells_in_subsampling = Donor_Original_Cells - Donor_Subsampled_Cells
        New_Cells_in_subsampling = Donor_Subsampled_Cells - Donor_Original_Cells
        
        New_Identities = Cell_identities_subsampled2.loc[Missing_Cells_in_subsampling]
        Becoming_Unassigned = New_Identities[New_Identities['donor_id']=='unassigned']
        
        New_Identities=New_Identities[New_Identities['donor_id']!='unassigned']
        Becoming_Doublet=New_Identities[New_Identities['donor_id']=='doublet']
        Becoming_Different_Donor=New_Identities[New_Identities['donor_id']!='doublet']
        
        # Now we add up all of this for the cells in the dataframe, that can be used for further quantification downstream.
        
        Reference_cell_identities2.loc[set(Becoming_Doublet.index),'Nr times becoming Doublet in subsampling'] = Reference_cell_identities2.loc[set(Becoming_Doublet.index),'Nr times becoming Doublet in subsampling']+1
        Reference_cell_identities2.loc[set(Becoming_Unassigned.index),'Nr times becoming Unassigned in subsampling'] = Reference_cell_identities2.loc[set(Becoming_Unassigned.index),'Nr times becoming Unassigned in subsampling']+1
        Reference_cell_identities2.loc[set(Becoming_Different_Donor.index),'Nr times becoming different donor in subsampling'] = Reference_cell_identities2.loc[set(Becoming_Different_Donor.index),'Nr times becoming different donor in subsampling']+1
        Reference_cell_identities2.loc[set(Becoming_Different_Donor.index),'New Donor Identities'] = Reference_cell_identities2.loc[Becoming_Different_Donor.index,'New Donor Identities']+Becoming_Different_Donor['donor_id']+';'
    Reference_cell_identities2['Nr vireo subsampling itterations']=Reference_cell_identities2['Nr vireo subsampling itterations']+1
Reference_cell_identities2.to_csv('subsampling_donor_swap_quantification.tsv',sep='\t')    

# Reference_cell_identities2[Reference_cell_identities2['Nr times becoming different donor in subsampling']==2]
print('Done')

