#!/usr/bin/env python
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-s", "--samplename", metavar="samplename", dest = "samplename",
                    help="", type=str,
                    default=None)

parser.add_argument("-v", "--variants_description", metavar="variants_description", dest = "variants_description",
                    help="", type=str,
                    default=None)

parser.add_argument("-s1", "--set1_uninformative_sites", metavar="set1_uninformative_sites", dest = "set1_uninformative_sites",
                    help="", type=str,
                    default=None)

parser.add_argument("-s2", "--set2_informative_sites", metavar="set2_informative_sites", dest = "set2_informative_sites",
                    help="", type=str,
                    default=None)

parser.add_argument("-p", "--positions_called_on", metavar="positions_called_on", dest = "positions_called_on",
                    help="", type=str,
                    default=None)

args = parser.parse_args()

positions_called_on = pd.read_csv(args.positions_called_on,sep='\t',header=None)
positions_called_on = set(positions_called_on[0].astype(str)+':'+positions_called_on[1].astype(str))
set2_informative_sites = pd.read_csv(args.set2_informative_sites,sep='\t',header=None)
set2_informative_sites = set(set2_informative_sites[0].astype(str)+':'+set2_informative_sites[1].astype(str))
set1_uninformative_sites = pd.read_csv(args.set1_uninformative_sites,sep='\t',header=None)
set1_uninformative_sites = set(set1_uninformative_sites[0].astype(str)+':'+set1_uninformative_sites[1].astype(str))
variants_description = pd.read_csv(args.variants_description,sep='\t')
samplename = args.samplename

variants_description['sample']= samplename
variants_description['total number of positions cellsp called on']= len(positions_called_on)
variants_description['informative positions cellsp called on']=len(positions_called_on.intersection(set2_informative_sites))
variants_description['uninformative positions cellsp called on']=len(positions_called_on.intersection(set1_uninformative_sites))
variants_description=variants_description.rename(columns={'total sites':'total sites attempted to call on','informative sites':'informative sites attempted to call on','uninformative sites':'uninformative sites attempted to call on'})
variants_description=variants_description[['sample', 'total sites attempted to call on','total number of positions cellsp called on','informative sites attempted to call on','informative positions cellsp called on','uninformative sites attempted to call on','uninformative positions cellsp called on']]
variants_description.to_csv(f"{samplename}_variants_description.tsv",sep='\t',index=False)
print('Done')