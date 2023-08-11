#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import glob

parser = argparse.ArgumentParser(
    description="""
        Combines all the Raw cellranger metadata to be passed in the h5ad
        """
)

parser.add_argument(
    '-d',
    action='store',
    dest='path',
    required=False,
    default='./results',
    help='Path')

options = parser.parse_args()
path = options.path

# If we run this code with pipeline that hasnt used genotypes we can not produce a cohort specific reports. therefore
if (os.path.exists(f"{path}/deconvolution/vireo_gt_fix")):
    try:
        try:
            prefix=f'{path}' #this is used for the updating reports posthoc, -> for this disable the next line
            project_name = os.listdir(f"{prefix}/handover/Summary_plots")[0]
            GT_MATCH = pd.read_csv(f"{path}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv",sep='\t')
            Donor_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Donor_Report.tsv",sep='\t')
            Tranch_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/{project_name}_Tranche_Report.tsv",sep='\t')
            try:
                Extra_Metadata_Donors = pd.read_csv(f"{prefix}/handover/Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
                Pipeline_inputs = pd.read_csv(f"{prefix}/Summary_plots/{project_name}/Fetch Pipeline/Input/input_table.tsv",sep='\t')
            except:
                Extra_Metadata_Donors = pd.read_csv(f"Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
                Pipeline_inputs = pd.read_csv(f"Summary_plots/{project_name}/Fetch Pipeline/Input/input_table.tsv",sep='\t')
        except:
            prefix='.'
            project_name = os.listdir(f"{prefix}/Summary_plots")[0]
            don_file = glob.glob(f"{path}/handover/Donor_Quantification_summary/*_Donor_Report.tsv")[0]
            tranch_file = glob.glob(f"{path}/handover/Donor_Quantification_summary/*_Tranche_Report.tsv")[0]
            GT_MATCH = pd.read_csv(f"{path}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv",sep='\t')
            Donor_Report = pd.read_csv(f"{don_file}",sep='\t')
            Tranch_Report = pd.read_csv(f"{tranch_file}",sep='\t')
            try:
                Extra_Metadata_Donors = pd.read_csv(f"{prefix}/Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
            except:
                try:
                    prefix=f'{path}/handover/'
                    Extra_Metadata_Donors = pd.read_csv(f"{prefix}/Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
                except:
                    _='No extra metadata provided'
                    prefix='.'
                    Extra_Metadata_Donors = pd.DataFrame(columns=['experiment_id'])
            Pipeline_inputs = pd.read_csv(f"{prefix}/Summary_plots/{project_name}/Fetch Pipeline/Input/input_table.tsv",sep='\t')
        # Split Donor Report by cohort


        Donor_Report2 = Donor_Report.set_index('Pool_ID.Donor_Id')
        
        try:
            Donor_Report2.insert(0,'Vacutainer ID','NONE')
        except:
            print('Vacutainer ID exists')
        try:
            Donor_Report2['Match Expected']='False'
        except:
            print('field doesnt exist')
        # here if we enable we can replace the expected number of samples if a mistake was made.
        # nr ukbb samples is already in there.

        # qc/Cardinal_45327_Jul_18_2022/work/62/aacd091ea7711fae264e884e122230/Summary_plots/Cardinal_45327_Jul_18_2022/Summary
        # t2 = GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline')]
        # GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline'),'Match Expected']=True
        GT_MATCH_CONFIDENT = GT_MATCH[GT_MATCH['Match Expected'] == True]

        def Generate_Report(GT_MATCH_CONFIDENT,pan):

            Total_Report = pd.DataFrame()               
            for i,row1 in GT_MATCH_CONFIDENT.iterrows():
                print(i)
                print(row1)
                print(f"{row1['pool']}_{row1['donor_query']}")
                
                Matched_Donor_report = ''
                Matched_Donor_report = Donor_Report[Donor_Report['Pool_ID.Donor_Id']==f"{row1['pool']}_{row1['donor_query']}"]
                
                Tranch_stats = Tranch_Report[Tranch_Report['Pool id'] ==row1['pool']]
                if ('celline' in row1['donor_gt']):
                    print('celline')
                    gt_match = row1['donor_gt'].split('_')[1].replace('celline_','')
                else:
                    gt_match = row1['donor_gt']
                
                Extra_Metadata_Donor1 = Extra_Metadata_Donors[Extra_Metadata_Donors['experiment_id'].str.contains(f"{row1['pool']}")]
                Extra_Metadata_Donor1 = Extra_Metadata_Donor1[Extra_Metadata_Donor1['experiment_id'].str.contains(f"{gt_match}")]
                if(len(Extra_Metadata_Donor1)>1):
                    print('More than 1 match')

                try:
                    Matched_Donor_report.insert(2, "Vacutainer ID", gt_match)
                except:
                    Matched_Donor_report['Vacutainer ID']=gt_match
                Matched_Donor_report['Experiment ID'] = project_name
                try:
                    Matched_Donor_report['Chromium channel number'] = Tranch_stats['Chromium channel number'].values[0]
                    Matched_Donor_report['Date of sample sequencing'] = Tranch_stats['Date of sample sequencing'].values[0]
                    cohort = row1['final_panel']
                    try:
                        live_cell_count = " or ".join(Extra_Metadata_Donor1['live_cell_count'])
                        viability =" or ".join(Extra_Metadata_Donor1['viability'])
                        SITE = " or ".join(Extra_Metadata_Donor1['SITE'])
                        Amount = " or ".join(Extra_Metadata_Donor1['Amount'])
                        RECIEVED = " or ".join(Extra_Metadata_Donor1['RECIEVED'])
                        try:
                            STATE = " or ".join(Extra_Metadata_Donor1['State'])
                        except:
                            STATE = 'NA'
                    except:
                        live_cell_count = 'NA'
                        viability ='NA'
                        SITE = 'NA'
                        Amount = 'NA'
                        RECIEVED = 'NA'
                        STATE = 'NA'
                        
                    try:
                        Matched_Donor_report.insert(8, 'lab_live_cell_count',live_cell_count)
                    except:
                        print('exists')
                        Matched_Donor_report['lab_live_cell_count'] =  live_cell_count

                    try:
                        Matched_Donor_report.insert(8, 'Sequencing time',STATE)
                    except:
                        print('exists')
                        Matched_Donor_report['Sequencing time'] =  STATE

                    try:
                        Matched_Donor_report.insert(8, 'viability', viability)
                    except:
                        print('exists')
                        Matched_Donor_report['viability'] =  viability

                    try:
                        Matched_Donor_report.insert(8, 'cohort', cohort)
                    except:
                        print('exists')
                        Matched_Donor_report['cohort'] =  cohort

                    try:
                        Matched_Donor_report.insert(8, 'Match Expected', row1['Match Expected'])
                    except:
                        print('exists')
                        Matched_Donor_report['Match Expected'] =  row1['Match Expected']

                    try:
                        Matched_Donor_report.insert(8, 'site', SITE)
                    except:
                        print('exists')
                        Matched_Donor_report['site'] =  SITE

                    try:
                        Matched_Donor_report.insert(8, 'amount recieved', Amount)
                    except:
                        print('exists')
                        Matched_Donor_report['amount recieved'] =  Amount

                    Matched_Donor_report['Date sample received'] =  RECIEVED
                    Matched_Donor_report['Donor id']=row1['donor_gt original'].split('_')[0]
                    Total_Report = pd.concat([Total_Report,Matched_Donor_report])
                except:
                    print('no successful matches in this cohort')
            return Total_Report

        GT_MATCH['final_panel']=GT_MATCH['final_panel'].str.replace("_[1-9][0-9]*$","", regex=True)
        for confident_panel in set(GT_MATCH['final_panel']):
            if (confident_panel!='NONE' and confident_panel!='GT_cell_lines'):
                print(confident_panel)
                pan=confident_panel.replace('GT_','')
                print(pan)
                GR_PANEL = GT_MATCH[GT_MATCH['final_panel'] == confident_panel]
                GT_CELLINE = GT_MATCH[GT_MATCH['final_panel'] == 'GT_cell_lines']
                for i,row1 in GT_CELLINE.iterrows():
                    celline = row1['donor_gt original'].split('_')[1]
                    if celline in row1['Good_ids expected']:
                        GT_CELLINE.loc[i,'Match Expected']=True

                GT_PANEL_CELLINE = pd.concat([GR_PANEL,GT_CELLINE])

                Total_Report = Generate_Report(GT_PANEL_CELLINE,pan)
                Total_Report2 = Total_Report.set_index('Pool_ID.Donor_Id')
                Total_Report = Total_Report.set_index('Pool_ID.Donor_Id')
                Donor_Report2.loc[Total_Report2.index,'Vacutainer ID']=Total_Report2.loc[Total_Report2.index,'Vacutainer ID']
                Donor_Report2.loc[Total_Report2.index,'Chromium channel number']=Total_Report2.loc[Total_Report2.index,'Chromium channel number']

                try: 
                    Donor_Report2.insert(12, "site",'NONE') 
                    Donor_Report2.insert(12, "lab_live_cell_count",'NONE') 
                    Donor_Report2.insert(12, "viability",'NONE') 
                    Donor_Report2.insert(5, "Match Expected",'False') 
                    Donor_Report2.insert(5, "Sequencing time",'False') 
                except: 
                    
                    print('exists')

                Donor_Report2.loc[Total_Report2.index,'site']=Total_Report2.loc[Total_Report2.index,'site']
                Donor_Report2.loc[Total_Report2.index,'lab_live_cell_count']=Total_Report2.loc[Total_Report2.index,'lab_live_cell_count']
                Donor_Report2.loc[Total_Report2.index,'viability']=Total_Report2.loc[Total_Report2.index,'viability']
                Donor_Report2.loc[Total_Report2.index,'Match Expected']=Total_Report2.loc[Total_Report2.index,'Match Expected']
                Donor_Report2.loc[Total_Report2.index,'Sequencing time']=Total_Report2.loc[Total_Report2.index,'Sequencing time']
                Donor_Report2 = Donor_Report2.reset_index().set_index('Pool ID')

                GT_MATCH2 = GT_MATCH.set_index('pool')
                GT_MATCH2=GT_MATCH2[['Emergency_ids expected','Good_ids expected']]
                GT_MATCH2 = GT_MATCH2.drop_duplicates().fillna("")
                GT_MATCH2['All IDs expected'] = GT_MATCH2['Emergency_ids expected']+','+GT_MATCH2['Good_ids expected']
                Donor_Report2['All IDs expected']=' '
                for idx1 in GT_MATCH2.index:
                    print(idx1)
                    Donor_Report2.loc[idx1,'All IDs expected']=GT_MATCH2.loc[idx1,'All IDs expected']
                Donor_Report2 = Donor_Report2.reset_index().set_index('Pool_ID.Donor_Id')

                if confident_panel == 'GT_ELGH':
                    confident_panel_name = 'Cardinal ELGH'
                    # confident_panel_name ='ELGH'
                elif confident_panel == 'GT_UKBB':
                    confident_panel_name = 'Cardinal UKB'

                Missing = pd.DataFrame(columns=['Pool'])
                Not_Expected = pd.DataFrame(columns=['Pool'])
                for id1 in set(GT_MATCH['pool']):
                    print(id1)
                    try:
                        Stats_File = pd.read_csv(f"{path}/gtmatch/{id1}/PiHAT_Stats_File_{id1}.csv",sep=',')
                        Stats_File['exp_id']=id1+'_'+Stats_File['donor_query']
                        Stats_File = Stats_File.set_index('exp_id')
                        overlapping_index=set(Total_Report.index).intersection(set(Stats_File.index))
                        if (len(overlapping_index)>0):
                            Total_Report.loc[overlapping_index,'PiHat: Expected']=Stats_File.loc[overlapping_index,'PiHat: Expected'].fillna('  ')
                            Total_Report.loc[overlapping_index,'Infered Relatednes (PiHAT>0.3)']=Stats_File.loc[overlapping_index,'Infered Relatednes (PiHAT>0.3)'].fillna('  ')
                            # Total_Report['Infered Relatednes (PiHAT>0.3)']=Stats_File['Infered Relatednes (PiHAT>0.3)'].fillna('  ')
                            # PiHat: Expected 
                            # if pan!='UKBB':
                            #     Total_Report.loc[overlapping_index,'PiHat: Expected']=Stats_File.loc[overlapping_index,'PiHat: Expected'].fillna('  ')
                            
                    except:
                        print('Stats file not available')
                    
                    Total_Report_samples = Total_Report[Total_Report['Pool ID']==id1]
                    Donor_Report2_samples = Pipeline_inputs[Pipeline_inputs['experiment_id']==id1]
                    Donor_Report2_samples['donor_vcf_ids']=Donor_Report2_samples['donor_vcf_ids'].fillna('')
                    DF1 = pd.DataFrame(Donor_Report2_samples.iloc[0]['donor_vcf_ids'].replace('\'','').split(','),columns=['col1'])
                    DF1['col1']= DF1['col1'].str.replace(r'U937_.*', 'U937', regex=True).astype('str')
                    DF1['col1']= DF1['col1'].str.replace(r'THP1_.*', 'THP1', regex=True).astype('str')
                    Expected_Samples = set(DF1.col1)

                    # print(Expected_Samples)
                    try:
                        Expected_Samples.remove('')
                    except:
                        _='cant'
                    try:
                        Expected_Samples2 = Extra_Metadata_Donors[Extra_Metadata_Donors.cohort == confident_panel_name]
                        Expected_Samples3 = Extra_Metadata_Donors[Extra_Metadata_Donors.cohort == 'Cardinal controls']
                        Expected_Samples3['donor']= Expected_Samples3['donor'].str.replace(r'U937_.*', 'U937', regex=True).astype('str')
                        Expected_Samples3['donor']= Expected_Samples3['donor'].str.replace(r'THP1_.*', 'THP1', regex=True).astype('str')
                        
                        combo = set(Expected_Samples2.donor).union(Expected_Samples3.donor)
                        Expected_Samples = set(combo).intersection(set(Expected_Samples))
                    except:
                        print('Donor level metadata with cohort info is not provided, assuming that all are there.')
                    All_Deconvouted_Samples = set(Total_Report_samples['Vacutainer ID'])
                    DF2= pd.DataFrame(Expected_Samples,columns=['col1'])
                    
                    # DF2['col1'] = DF2['col1'].str.replace(r'U937_.*', 'U937', regex=True).astype('str')
                    # DF2['col1'] = DF2['col1'].str.replace(r'THP1_.*', 'THP1', regex=True).astype('str')
                    DF3= pd.DataFrame(All_Deconvouted_Samples,columns=['col1'])
                    
                    Expected_Samples = set(DF2.col1.str.replace('^0*', ''))
                    All_Deconvouted_Samples = set(DF3.col1.str.replace('^0*', ''))
                    Missing_Samples = set(Expected_Samples)-set(All_Deconvouted_Samples)
                    Not_Expected_Samples = set(All_Deconvouted_Samples)-set(Expected_Samples)
                    if (len(Missing_Samples)>0):
                        Missing_Samples2=pd.DataFrame(Missing_Samples,columns=['Sample'])
                        Missing_Samples2['Pool']=id1
                        Missing=pd.concat([Missing,Missing_Samples2])

                        Not_Expected_Samples2=pd.DataFrame(Not_Expected_Samples,columns=['Sample'])
                        Not_Expected_Samples2['Pool']=id1
                        Not_Expected=pd.concat([Not_Expected,Not_Expected_Samples2])
                    else:
                        print('all good')
                Not_Expected = Not_Expected.set_index('Pool')
                Missing = Missing.set_index('Pool')
                try:
                    Donor_Report2.loc[Total_Report.index,'PiHat: Expected']=Total_Report.loc[Total_Report.index,'PiHat: Expected']
                    Donor_Report2.loc[Total_Report.index,'Infered Relatednes (PiHAT>0.3)']=Total_Report.loc[Total_Report.index,'Infered Relatednes (PiHAT>0.3)']
                except:
                    _=''

                if (pan=='UKBB'):
                    # For UKB samples we only return the expected samples
                    Total_Report = Total_Report[Total_Report['Match Expected']]
                else:
                    # Disable this if we dont want to limit to only expected
                    Total_Report = Total_Report[Total_Report['Match Expected']]



                
                
                # Missing_Samples = set(Expected_Samples.donor)-set(Total_Report['Vacutainer ID'])
                # This should loop through each of the inputs - expected ws deconvoluted.
                try:
                    os.mkdir(f'Summary_plots/{project_name}/Summary/{pan}_REPORT')
                except:
                    print('Dir exists')
                # if (len(Missing)>0):
                Missing.to_csv(f'Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_Missing_{pan}_Donors.tsv',sep='\t')
                # if (len(Not_Expected)>0):
                Not_Expected.to_csv(f'Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_Not_Expected_{pan}_Donors.tsv',sep='\t')

                Total_Report.to_csv(f'Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_{pan}_Report.tsv',sep='\t')
            
        Donor_Report2.to_csv(f"results/handover/Donor_Quantification_summary/{project_name}_Donor_Report.tsv",sep='\t',index=True)
        Donor_Report2.to_csv(f"Summary_plots/{project_name}/Summary/{project_name}_Donor_Report.tsv",sep='\t',index=True)

        print('Done')
    except:
        print('Pipeline most likely run without genotypes, therefore no cohorts available, skipping processing')
else:
    print('Pipeline run without genotypes, therefore no cohorts available, skipping processing')