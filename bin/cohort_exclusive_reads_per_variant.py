#!/usr/bin/env python3

__date__ = '2023-07-13'
__version__ = '0.0.1'
import argparse
import pickle 
import pandas as pd
import random
import numpy as np
import multiprocessing as mp
from multiprocessing import Lock
import os
import gzip
import time

class VCF_Loader:
    
    def __init__(self, vcf_file, biallelic_only=True,
                        sparse=False, format_list=['GT']):
        self.vcf_file = vcf_file
        self.load_sample = True
        self.biallelic_only = biallelic_only
        self.sparse = sparse
        self.record_dict={}
        self.reset()
        self.format_list = format_list
        self.exclusive_donor_variants = {}
        self.curently_pushing =[] #this is a lock value to check if rhe curent field is updated so to avaid the race for update
        self.last_count=-1
        self.reset_c()
    
    def reset_c(self):
        self.record_times=0
        
    def reset(self):
        self.exclusive_donor_variants ={}
                
    def myfunc(self):
        print(f"Hello my name is {self.biallelic_only}" )
        
    def load_sample_mp(self,line,obs_ids,count,format_list):
        '''
        takes VCF lines and extracts all format fields for those where GT !='.'
        '''
        list_val = line.rstrip().split("\t") #[:5] #:8
        idx = find(list_val[8].split(':'),'GT')[0]#find index of GT field as GT will tell us what variants are called
        if len(list_val[3]) > 1 or len(list_val[4]) > 1:
            # CURRENTLY DEALS ONLY WITH BIALELIC
            print(f'{idx} var not bialelic')
        else:
            list_val2 = list_val[9:]
            obs = pd.DataFrame(obs_ids)
            lv = pd.DataFrame(list_val2)
            lv_proc =lv[0].str.split(':').str[idx]
            gt_exists = lv_proc[lv_proc != '.']
            idx2 = gt_exists.index
            obs_with_gt = obs.loc[idx2.values]
            obs_with_gt = list(obs_with_gt[0].values)
            list_val_with_gt = lv.loc[idx2.values]
            list_val_with_gt = list(list_val_with_gt[0].values)
            random.seed(count)
            c = list(zip(obs_with_gt, list_val_with_gt))
            random.shuffle(c)
            obs_with_gt, list_val_with_gt = zip(*c)
            # self.append_results([obs_with_gt,list_val_with_gt,idx,list_val,count])

        return [obs_with_gt,list_val_with_gt,idx,list_val,count,format_list]#add format_list to the return value as we need this for the next step


    def set_results(self,to_set,id):
        # Recod to disk to save the loading mmeory time.
        with open(f'tmp_{id}.pkl', 'wb') as f:
            pickle.dump(to_set, f)
        self.record_dict[id]=f'tmp_{id}.pkl'
    

    def append_results(self,result):
        # exclusive_donor_variants
        obs_with_gt= result[0]
        list_val_with_gt= result[1]
        idx = result[2]
        list_val = result[3]
        count = result[4]
        format_list = result[5]#list of required format fields
        #get indexes of required format fields (apart from GT which has already been taken care of)
        additional_field_idxs = []
        for fmt in format_list:
            if not fmt == 'GT':
                idx_addn = find(list_val[8].split(':'), fmt)[0]
                additional_field_idxs.append(idx_addn)
        # print(additional_field_idxs)
        # exit(0)

        count11=0
        # r = random.random()
        # Issue is that this slows down after number of entries is recorded. So recoding takes longer and longer.
        # every 500 itterations we push the data to a dictionary, later we combine these together.
        if (count % 200 == 0):
            print(f'recording and resetting memory {count}')
            # self.record_dict[count]=self.exclusive_donor_variants
            self.set_results(self.exclusive_donor_variants,count)
            self.reset()  
            self.reset_c()        
        
        for ob_id in obs_with_gt:
            donor_loc_in_list = count11
            alleles = list_val_with_gt[donor_loc_in_list].split(':')[idx]
            #append any additional format fields to alleles
            if len(additional_field_idxs) > 0:
                for idx_addnl in additional_field_idxs:
                    fmt_val = list_val_with_gt[donor_loc_in_list].split(':')[idx_addnl]
                    alleles = alleles + '_' + fmt_val

            if not alleles.startswith('.'):
                ids = "_".join([list_val[x] for x in [0, 1, 3, 4]])
                donor_var = f"{ids}_{alleles}"
                while ob_id in self.curently_pushing:
                    time.sleep(r*0.01)
                self.curently_pushing.append(ob_id)           
                try:
                    self.exclusive_donor_variants[ob_id].add(donor_var)
                    self.record_times=self.record_times+1
                except:
                    self.exclusive_donor_variants[ob_id]=set()
                    self.exclusive_donor_variants[ob_id].add(donor_var)
                    self.record_times=self.record_times+1
                self.curently_pushing.remove(ob_id)
                # self.exclusive_donor_variants['CTGAAACGTAAGTTCC-1']
            count11+=1 

    def combine_written_files(self):#this is for VCF loader class
        to_export = self.exclusive_donor_variants
        for val1 in self.record_dict.values():
            # here remove the int files.
            print(f"merging temp file: {val1}")
            with open(val1, 'rb') as f:
                loaded_dict = pickle.load(f)
                for k1 in loaded_dict.keys():
                    try:
                        to_export[k1]=to_export[k1].union(loaded_dict[k1])
                    except:
                        to_export[k1]=set()
                        to_export[k1]=to_export[k1].union(loaded_dict[k1])
            os.remove(val1)
        return to_export
    
    
    def load_VCF_batch_paralel(self):
        """
        Load whole VCF file by utilising multiple cores to speed up loading of large cell files
        -------------------
        Initially designed to load VCF from cellSNP output, requiring 
        1) all variants have the same format list;
        2) a line starting with "#CHROM", with sample ids.
        If these two requirements are satisfied, this function also supports general
        VCF files, e.g., genotype for multiple samples.

        Note, it may take a large memory, please filter the VCF with bcftools first.
        """
        
        vcf_file = self.vcf_file
        biallelic_only = self.biallelic_only
        load_sample= self.load_sample
        sparse = self.sparse
        format_list= self.format_list
        pool = mp.Pool(cpus)
        
        
        import time
        if vcf_file[-3:] == ".gz" or vcf_file[-4:] == ".bgz":
            infile = gzip.open(vcf_file, "rb")
            is_gzip = True
        else:
            infile = open(vcf_file, "r")
            is_gzip = False
        
        FixedINFO = {}
        contig_lines = []
        comment_lines = []
        var_ids, obs_ids, obs_dat = [], [], []
        count=0 #57077    
        for line in infile:
            count+=1
            # if count>10000:
            #     break
            if is_gzip:
                line = line.decode('utf-8')
            if line.startswith("#"):
                if line.startswith("##contig="):
                    contig_lines.append(line.rstrip())
                if line.startswith("#CHROM"):
                    if load_sample:
                        obs_ids = line.rstrip().split("\t")[9:]
                        for ob_id in obs_ids:
                            self.exclusive_donor_variants[ob_id]=set()
                    key_ids = line[1:].rstrip().split("\t")[:8]
                    for _key in key_ids:
                        FixedINFO[_key] = []
                else:
                    comment_lines.append(line.rstrip())
            else:
                pool.apply_async(self.load_sample_mp, args=([line,obs_ids,count,format_list]),callback=self.append_results)
                del line
        self.last_count=count
        pool.close()
        pool.join()
        
        output = self.combine_written_files()
        return output
    

class Read_Counts_Per_Var:
    def __init__(self, cell_variants_ad, cell_assignments_table, donor_assignments_table, private_vars):
        #self.reset()
        self.donor_assignments_table=donor_assignments_table
        self.cell_assignments_table=cell_assignments_table
        self.cell_variants=cell_variants_ad
        self.private_vars = private_vars
        self.result_dict={}


    # def reset(self):
    #     self.read_counts ={}


    # def set_results(self,to_set,id):
    #     # Recod to disk to save the loading mmeory time.
    #     with open(f'tmp_{id}.pkl', 'wb') as f:
    #         pickle.dump(to_set, f)
    #     self.record_dict[id]=f'tmp_{id}.pkl'


    # def combine_written_files(self):#this one is for concordance class
    #     to_export = self.read_counts
    #     for val1 in self.record_dict.values():
    #         # here remove the int files.
    #         print(f"merging temp file: {val1}")
    #         with open(val1, 'rb') as f:
    #             loaded_dict = pickle.load(f)
    #             for k1 in loaded_dict.keys():
    #                 to_export[k1]=loaded_dict[k1]
    #         os.remove(val1)
    #     return to_export
    

    def vars_to_pandas(self,vars):
        vars_df = pd.DataFrame(vars)
        if len(vars_df) > 0:
            split_str=vars_df[0].str.split("_")
            vars_df['ids'] = split_str.str[0]+'_'+split_str.str[1]+'_'+split_str.str[2]+'_'+split_str.str[3]
            vars_df['pos'] = split_str.str[0]+'_'+split_str.str[1]
            vars_df['vars'] = split_str.str[4]
            vars_df['vars'] = vars_df['vars'].str.replace('|','/',regex=False)
            vars_df = vars_df[vars_df['vars']!='./.']
            vars_df.loc[vars_df['vars']=='0/1','vars']='1/0'
            vars_df['combo']= vars_df['ids']+'_'+vars_df['vars']
        return vars_df
    

    def read_counts_per_var_table_production(self, cell_vars, cohort):
        """
        Identify UKB/ELGH specific vars in cell
        """
        private_vars = self.private_vars
        cell_vars_df = self.vars_to_pandas(cell_vars)
        cell_vars_df['AD'] = cell_vars_df[0].str.split("_").str[5].astype(int)
        cell_vars_df['DP'] = cell_vars_df[0].str.split("_").str[6].astype(int)
        cell_vars_df['OTH'] = cell_vars_df[0].str.split("_").str[7].astype(int)
        #sum reads in cell vars
        cell_vars_df['total_dp'] = cell_vars_df['DP'] + cell_vars_df['OTH']
        #intersect cell vars with private vars
        cell_private_vars = cell_vars_df[cell_vars_df['ids'].isin(private_vars)]
        # print(self.result_dict)
        # print(cell_private_vars)
        #now iterate through cell_private_vars by row and add the AD/total to the approriate dict items
        for index, row in cell_private_vars.iterrows():
            idr = row['ids']
            adr = row['AD']
            totr = row['total_dp']
            if cohort == 'ELGH':
                self.result_dict[idr]['AD_ELGH'] = self.result_dict[idr]['AD_ELGH'] + adr
                self.result_dict[idr]['total_reads_ELGH'] = self.result_dict[idr]['total_reads_ELGH'] + totr
            elif cohort == 'UKB':
                self.result_dict[idr]['AD_UKB'] = self.result_dict[idr]['AD_UKB'] + adr
                self.result_dict[idr]['total_reads_UKB'] = self.result_dict[idr]['total_reads_UKB'] + totr                
        #     print(self.result_dict[idr])
        # exit(0)

    # def append_results_read_counts(self, result):
    #     count=result[13]
    #     print(count)

    #     # self.read_counts[f'{result[0]} --- {result[1]}'] = {'cell':result[0],
    #     #                                                     'donor':result[1],
    #     #                                                     'elgh_ad':result[2],
    #     #                                                     'ukb_ad':result[3],
    #     #                                                     'elgh_var_total_reads': result[4],
    #     #                                                     'elgh_var_total_reads_no_hom_ref': result[5],
    #     #                                                     'ukb_var_total_reads': result[6],
    #     #                                                     'ukb_var_total_reads_no_hom_ref': result[7],
    #     #                                                     'elgh_median_vaf': result[8],
    #     #                                                     'ukb_median_vaf': result[9],
    #     #                                                     'cellsnp_read_count': result[10],
    #     #                                                     'cellsnp_median_vaf': result[11],
    #     #                                                     'cohort':result[12]}

    #     if (count % 200 == 0):
    #         print(f'recording and resetting memory {count}')
    #         # self.record_dict[count]=self.exclusive_donor_variants
    #         self.set_results(self.read_counts,count)
    #         self.reset()  
    #     _=""

    
    # def set_results(self,to_set,id):
    #     # Recod to disk to save the loading mmeory time.
    #     with open(f'tmp_{id}.pkl', 'wb') as f:
    #         pickle.dump(to_set, f)
    #     self.record_dict[id]=f'tmp_{id}.pkl'


    def results_table(self):
        donor_assignments_table=self.donor_assignments_table
        cell_assignments_table=self.cell_assignments_table
        cell_variants_ad=self.cell_variants
        private_vars = self.private_vars

        #set up results dict
        for v in private_vars:
            self.result_dict[v] = {'AD_ELGH': 0, 'total_reads_ELGH': 0, 'AD_UKB': 0, 'total_reads_UKB': 0}
        
        # pool = mp.Pool(cpus)
        # count = 0
        for i,row1 in donor_assignments_table.iterrows():
            donor_in_question = row1['donor_query']
            donor_gt_match = row1['donor_gt']

            if (donor_gt_match=='NONE'):
                continue

            cohort = 'UNKNOWN'
            donor_split = donor_gt_match.split("_")
            if (len(donor_split) == 2) and (donor_split[0] == donor_split[1]):
                cohort = 'UKB'
            elif (len(donor_split) == 3) and (len(donor_split[0]) == 14):
                cohort = 'ELGH'
            if cohort == 'UNKNOWN':#we don't want any cells that are not ELGH or UKB
                continue

            Cells_to_keep_pre = list(set(cell_assignments_table.loc[cell_assignments_table['donor_id']==donor_in_question,'cell']))
                
            for cell1 in Cells_to_keep_pre:
                #count+=1
                cell_vars = cell_variants_ad[cell1]#select variants for just one cell
                #self.read_counts[f'{cell1} --- {donor_gt_match}']={}
                # print(cell_vars)
                # exit(0)
                self.read_counts_per_var_table_production(cell_vars, cohort)
                # self.append_results_read_counts(result)

        return self.result_dict
        # exit(0)
        # pool.close()
        # pool.join()
        # output = self.combine_written_files()
        # return output

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpus', action='store', required=True, type=int)
    parser.add_argument('--cell_vcf', action='store', required=True)
    parser.add_argument('--cell_assignments', action='store', required=True)
    parser.add_argument('--donor_assignments', action='store', required=True)
    parser.add_argument('--ukb_only_vars', action='store', required=True)
    parser.add_argument('--elgh_only_vars', action='store', required=True)
    parser.add_argument('--outfile_prefix', action='store', required=True)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    return args


def find(lst, a):
    return [i for i, x in enumerate(lst) if x==a ]


def get_vars_from_file(vars_file):
    """
    create a set from a file
    """
    sites = set()
    with open(vars_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            var = l.rstrip()
            sites.add(var)
    return sites


if __name__ == "__main__":
    options = get_options()
    cpus = options.cpus
    cell_vcf=options.cell_vcf
    cell_assignments=options.cell_assignments
    donor_assignments = options.donor_assignments
    ukb_only_vars_file = options.ukb_only_vars
    elgh_only_vars_file = options.elgh_only_vars
    outfile_prefix = options.outfile_prefix

    elgh_outfile = outfile_prefix + "_ELGH_private_vars.txt"
    ukb_outfile = outfile_prefix + "_UKB_private_vars.txt"

    ukb_only_vars = get_vars_from_file(ukb_only_vars_file)
    elgh_only_vars = get_vars_from_file(elgh_only_vars_file)

    donor_assignments_table = pd.read_csv(donor_assignments)
    cell_assignments_table = pd.read_csv(cell_assignments,sep='\t')

    if options.debug:
        with open('tmp_cell_variants.pkl', 'rb') as f:
            cell_variants = pickle.load(f) 
    else:
        print('---Loading cell VCF----')
        loader1 = VCF_Loader(cell_vcf, biallelic_only=True,
                            sparse=False, format_list=['GT', 'AD', 'DP', 'OTH'])
        cell_variants = loader1.load_VCF_batch_paralel()
        del loader1

        with open(f'tmp_cell_variants.pkl', 'wb') as f:
            pickle.dump(cell_variants, f)


    ukb_read_counts = Read_Counts_Per_Var(cell_variants, cell_assignments_table, donor_assignments_table, ukb_only_vars)
    output_table_ukb = ukb_read_counts.results_table()

    result_ukb = pd.DataFrame.from_records(output_table_ukb).T
    result_ukb.to_csv(ukb_outfile,sep='\t')

    elgh_read_counts = Read_Counts_Per_Var(cell_variants, cell_assignments_table, donor_assignments_table, elgh_only_vars)
    output_table_elgh = elgh_read_counts.results_table()

    result_elgh = pd.DataFrame.from_records(output_table_elgh).T
    result_elgh.to_csv(elgh_outfile,sep='\t')
    print('Processing Done')