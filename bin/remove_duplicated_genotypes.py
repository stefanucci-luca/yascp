#!/usr/bin/env python

__date__ = '2023-06-26'
__version__ = '0.0.2'

import pandas as pd
import argparse
import pandas as pd
import re
# Vireo is not working well if multiple same donor genotypes are provided.
# For this reasone we have to ensure they are removed.
parser = argparse.ArgumentParser(
    description="""
        Joins all the outputs from the GT match in an expected input format for subsequent relatedness estimation.
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-d', '--donors',
    action='store',
    dest='donors',
    required=True,
    help='donors'
)

parser.add_argument(
    '-b', '--bridge',
    action='store',
    dest='bridge',
    required=True,
    help='bridge'
)

parser.add_argument(
    '-i', '--yascp_input_data',
    action='store',
    dest='yascp_input_data',
    required=True,
    help='yascp_input_data'
)

parser.add_argument(
    '-n', '--name',
    action='store',
    dest='name',
    required=True,
    help='name'
)

options = parser.parse_args()
donors = pd.read_csv(options.donors,sep='\t',header=None)
bridge = pd.read_csv(options.bridge,sep='\t')
input_table=pd.read_csv(options.yascp_input_data,sep='\t')
pool_name = options.name
input_table=pd.read_csv(options.yascp_input_data,sep='\t')
input_expected = input_table.loc[input_table['experiment_id']==pool_name,'donor_vcf_ids'].values[0]
# /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Cardinal_46499_Jan_21_2023/results/yascp_inputs/input.tsv/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Cardinal_46499_Jan_21_2023/results/yascp_inputs/input.tsv
newlist = [] # empty list to hold unique elements from the list
duplist = []
mylist = bridge['oragene_id']
for i in mylist:
    if i not in newlist:
        newlist.append(i)
    else:
        duplist.append(i) 

donors['mapping']=''
for don1 in donors[0]:
    # print(don1)
    mappings = bridge[bridge['oragene_id']==don1]['s00046_id']
    
    if len(mappings)>0:
        mapping=mappings.values[0]
        for mapping in mappings:
            if mapping in input_expected:
                donors.loc[donors[0]==don1,'mapping']=mapping
                # there are cases where two genotype ids map to same expected donor. In this case we should pick the one expected.
                continue
    mappings = bridge[bridge['oragene_id']==don1.split('_')[0]]['s00046_id']
    if len(mappings)>0:
        mapping=mappings.values[0]
        for mapping in mappings:
            if mapping in input_expected:
                donors.loc[donors[0]==don1,'mapping']=mapping
                # there are cases where two genotype ids map to same expected donor. In this case we should pick the one expected.
                continue
        
donors = donors.sort_values('mapping')
donors = donors.drop_duplicates(subset=['mapping'])
donors[0].to_csv('t.tsv',index=False,header=None)
expected = set(input_expected.replace("'",'').split(','))
missing = expected-set(donors['mapping'])
print('Done')

