#!/usr/bin/env python3

__date__ = '2023-04-14'
__version__ = '0.0.1'

import argparse
import random
import pickle 
import pandas as pd
import gzip
import numpy as np
import time
import multiprocessing as mp
from multiprocessing import Lock
import os

use_only_informative_snps=True #When this is set to True we ignore any sites that have no difference between all the individuals.
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
    
    def find(self,lst, a):
        return [i for i, x in enumerate(lst) if x==a ]
    
    def reset_c(self):
        self.record_times=0
        
    def reset(self):
        self.exclusive_donor_variants ={}
                
    def myfunc(self):
        print(f"Hello my name is {self.biallelic_only}" )
        
    def load_sample_mp(self,line,obs_ids,count,format_list):
        list_val = line.rstrip().split("\t") #[:5] #:8
        idx = self.find(list_val[8].split(':'),format_list[0])[0]
        coment_fields = list_val[7].split(';')
        
        if len(list_val[3]) > 1 or len(list_val[4]) > 1:
            # CURRENTLY DEALS ONLY WITH BIALELIC
            _=f'{idx} var not bialelic'
            # print(f'{idx} var not bialelic')
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
            test_if_there_is_difference = set(pd.DataFrame(list_val_with_gt)[0].str.split(":").str[0])
            nr_genotypes = len(test_if_there_is_difference)
            if use_only_informative_snps==True:
                if nr_genotypes>1:
                    return [obs_with_gt,list_val_with_gt,idx,list_val,count,coment_fields]
                else:
                    _='Site is not informative'
                    # print('Site is not informative')
            else:
                return [obs_with_gt,list_val_with_gt,idx,list_val,count,coment_fields]
        

    def set_results(self,to_set,id):
        # Recod to disk to save the loading mmeory time.
        with open(f'tmp_{id}.pkl', 'wb') as f:
            pickle.dump(to_set, f)
        self.record_dict[id]=f'tmp_{id}.pkl'
    
    def append_results(self,result):
        # exclusive_donor_variants
        if result!=None:
            obs_with_gt= result[0]
            list_val_with_gt= result[1]
            idx = result[2]
            list_val = result[3]
            count = result[4]
            comment_fields= result[5]
            # if any("R2=" in s for s in comment_fields):
            #     comment_fields
            # print(count)
            # print(self.record_times)
            # if
            # print(self.last_count)
            count11=0
            # r = random.random()
            # Issue is that this slows down after number of entries is recorded. So recoding takes longer and longer.
            # every 500 itterations we push the data to a dictionary, later we combine these together.
            if (count % 300 == 0):
                print(f'recording and resetting memory {count}')
                # self.record_dict[count]=self.exclusive_donor_variants
                self.set_results(self.exclusive_donor_variants,count)
                self.reset()  
                self.reset_c()        
            
            for ob_id in obs_with_gt:
                donor_loc_in_list = count11
                alleles = list_val_with_gt[donor_loc_in_list].split(':')[idx]
                if alleles!='.':
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

    def combine_written_files(self):
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
                if cpus==1:
                    resu = self.load_sample_mp(line,obs_ids,count,format_list)
                    self.append_results(resu)
                else:
                    pool.apply_async(self.load_sample_mp, args=([line,obs_ids,count,format_list]),callback=self.append_results)
                del line
        self.last_count=count
        pool.close()
        pool.join()
        
        output = self.combine_written_files()
        
        return output

def norm_genotypes(expected_vars):
    expected_vars = pd.DataFrame(expected_vars)
    expected_vars['ids'] = expected_vars[0].str.split("_").str[:-1].str.join('_')
    expected_vars['pos'] = expected_vars[0].str.split("_").str[:2].str.join('_')
    expected_vars['vars'] = expected_vars[0].str.split("_").str[-1].str.join('')
    expected_vars['vars'] = expected_vars['vars'].str.replace('|','/',regex=False)
    expected_vars = expected_vars[expected_vars['vars']!='./.']
    expected_vars.loc[expected_vars['vars']=='0/1','vars']='1/0'
    expected_vars['combo']= expected_vars['ids']+'_'+expected_vars['vars']
    return expected_vars

def donor_exclusive_sites(exclusive_don_variants2):
    # Here we generate a function for determining the sites that are donor exclusive
    donor_distinct_sites = {}
    for col1 in exclusive_don_variants2.keys():
        comparisons =[]
        to_compare = []
        for col2 in exclusive_don_variants2.keys():
            if col1==col2:
                # we set this as the unique entry
                # print('1')
                to_compare = set(exclusive_don_variants2[col2])
            else:
                # We combine all the variants in one list
                comparisons+=list(exclusive_don_variants2[col2])
                # print('2')
        # print('comparison')  
        comparisons_all = set(comparisons) 
        comparisons_all_norm = norm_genotypes(comparisons_all)
        comparisons_all=set(comparisons_all_norm['combo'])
        
        to_compare = set(to_compare)
        to_compare_norm = norm_genotypes(to_compare)
        to_compare=set(to_compare_norm['combo'])
        # Make sure we account for hap types - phased/unphased 
        distinct_donor_sites = to_compare - comparisons_all
        donor_distinct_sites[col1]=distinct_donor_sites
        # Perform the distinct set function.
    return donor_distinct_sites

if __name__ == "__main__":
    # cpus=6
    # vcf_file='/lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra1/work/cf/5788c53e50951eb949dcb8f530be48/Study_Merge_AllExpectedGT_5PTFD3NBJ_out.vcf.gz'
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Runs BBKNN. Assumes PCS have already been
            calculated.
            """
    )
    parser.add_argument(
        '-cpus', '--cpus',
        action='store',
        dest='cpus',
        required=True,
        help='cpus'
    ) 
    
    parser.add_argument(
        '-vcf', '--vcf',
        action='store',
        dest='vcf',
        required=True,
        help='vcf'
    )     

    parser.add_argument(
        '-cellsnp', '--cellsnp',
        action='store',
        dest='cellsnp',
        required=True,
        help='cellsnp'
    )   

    options = parser.parse_args()
    vcf = options.vcf
    cellsnp = options.cellsnp
    cpus = int(options.cpus)

    print('---Genotype loader init----')    
    loader2 = VCF_Loader(vcf, biallelic_only=True,
                    sparse=False, format_list=['GT'])

    GT_Matched_variants = loader2.load_VCF_batch_paralel()
    # donor_distinct_sites = donor_exclusive_sites(GT_Matched_variants)
    cellsnp = pd.read_csv(cellsnp,sep='\t',comment='#',header=None)
    cellsnp[0]=cellsnp[0].astype(str)
    cellsnp[1]=cellsnp[1].astype(str)
    donor_informaive_sites = GT_Matched_variants
    ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    # chr10	20733	.	T	A	.	.	AF_UKBB=.;AF_ELGH=.;AF_spike_in=0.594649
    # chr10	44822	.	A	C	.	.	AF_UKBB=0.0452576;AF_ELGH=.;AF_spike_in=0.159744
    # chr10	45402	.	T	C	.	.	AF_UKBB=.;AF_ELGH=0.654687;AF_spike_in=.
    exta_snps = pd.DataFrame()
    # in the folowing
    iteration_dataframe = donor_informaive_sites #Here we provide either donor distinct sites or donor informative sites. 
    for key1 in iteration_dataframe.keys():
        all_sites = pd.DataFrame(iteration_dataframe[key1],columns=['f1'])
        splits = all_sites['f1'].str.split('_')
        all_sites['#CHROM']=splits.str[0]
        all_sites['POS']=splits.str[1]
        all_sites['ID']=f'.'
        all_sites['REF']=splits.str[2]
        all_sites['ALT']=splits.str[3]
        all_sites['QUAL']=f'.'
        all_sites['FILTER']=f'.'
        all_sites['INFO']=f'.'
        all_sites = all_sites.drop(columns=['f1'])
        exta_snps=pd.concat([exta_snps,all_sites])
    # All_Extra_informative_Sites = exta_snps.drop_duplicates()
    exta_snps.columns = cellsnp.columns   
    exta_snps=exta_snps.drop_duplicates(subset=[0, 1])
    cellsnp_exta_snps=pd.concat([cellsnp,exta_snps])
    exta_snps.index=exta_snps[0]+'_'+exta_snps[1]
    cellsnp.index=cellsnp[0]+'_'+cellsnp[1]
    set1_uninformative_sites=set(cellsnp.index)-set(exta_snps.index)
    informative_sites_covered_in_default_panel = set(exta_snps.index)-set(cellsnp.index)
    
    cellsnp_exta_snps = cellsnp_exta_snps.drop_duplicates(subset=[0, 1])
    set1_uninformative_sites = cellsnp.loc[set1_uninformative_sites]
    set2_informative_sites =exta_snps
    description = pd.DataFrame([{'total sites':len(cellsnp_exta_snps),'informative sites':len(set2_informative_sites),'uninformative sites':len(set1_uninformative_sites),'informative sites covered in initial panel':len(informative_sites_covered_in_default_panel)}])
    description.to_csv('variants_description.tsv',sep='\t',index=False)
    cellsnp_exta_snps.to_csv('cellsnp_variants.tsv',sep='\t',index=False,header=False)
    set1_uninformative_sites.to_csv('set1_uninformative_sites.tsv',sep='\t',index=False,header=False)
    set2_informative_sites.to_csv('set2_informative_sites.tsv',sep='\t',index=False,header=False)
    print('Done')
    
    


