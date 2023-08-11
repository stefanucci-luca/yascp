#!/usr/bin/env python3
import pandas as pd
import argparse
__date__ = '2021-11-04'
__version__ = '0.0.1'

parser = argparse.ArgumentParser(
    description="""
        Subsetting to 80% variants
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-random_state', '--random_state',
    action='store',
    dest='random_state',
    default=1,type=int
)

parser.add_argument(
    '-vcf', '--vcf',
    action='store',
    dest='vcf',
    default=None,
)

parser.add_argument(
    '-r', '--rate',
    action='store',
    dest='rate',
    default=80,type=int
)

parser.add_argument(
    '-o', '--order',
    action='store',
    dest='order',
    default=None,type=str
)


options = parser.parse_args()
rate = options.rate/100


VCF = pd.read_csv(options.vcf,sep='\t',comment='#',header=None)
VCF2 = VCF.sample(frac=rate,random_state=options.random_state)
VCF2=VCF2.sort_index()
VCF2['match_name'] = VCF2[0].astype(str)+';'+VCF2[1].astype(str)+';'+VCF2[2].astype(str)+';'+VCF2[3].astype(str)+';'+VCF2[4].astype(str)+';'+VCF2[5].astype(str)+';'+VCF2[6].astype(str)+';'+VCF2[7].astype(str)

# Here we make sure that the order is the same as in the base csv, otherwise sample assignments are wrong in deconvolutions.
vcf_order = pd.read_csv(options.order,sep='\t',comment='#',header=None)
vcf_order['order_index'] = vcf_order.index
vcf_order['match_name'] = vcf_order[0].astype(str)+';'+vcf_order[1].astype(str)+';'+vcf_order[2].astype(str)+';'+vcf_order[3].astype(str)+';'+vcf_order[4].astype(str)+';'+vcf_order[5].astype(str)+';'+vcf_order[6].astype(str)+';'+vcf_order[7].astype(str)
vcf_order=vcf_order.set_index('match_name')
VCF2 = VCF2.set_index('match_name')
# Now sort based on order index and emit the corect matrix.
VCF2['order_index']= vcf_order['order_index']
VCF2 = VCF2.sort_values(by=['order_index'])
VCF2 = VCF2.drop('order_index',axis=1)




VCF2.to_csv('random_variants.tsv',sep='\t',index=False,header=False)
print('Done')