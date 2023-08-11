#!/usr/bin/env python

__date__ = '2021-01-20'
__version__ = '0.0.1'

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# This code combines datasets in one file and plots cell concordances and discordances.
# There will always be less cells in the cell_concordance table since there are cells that do not have a donor assignments:
# the ones that have resulted in a gt check of:
# NONE
# Unassigned
# Doublet


"""Run CLI."""
parser = argparse.ArgumentParser(
    description="""
        Combines data and plots per pool statistics
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-cc', '--Cell_Concordance',
    action='store',
    dest='cc',
    required=True,
    help='Cell_Concordance'
)

parser.add_argument(
    '-sq', '--Swap_Quant',
    action='store',
    dest='sq',
    required=True,
    help='Swap_Quant'
)

parser.add_argument(
    '-name', '--name',
    action='store',
    dest='name',
    required=True,
    help='name.'
)

parser.add_argument(
    '-run', '--run',
    action='store',
    dest='run',
    required=False,
    default=None,
    help='run name'
)

options = parser.parse_args()

cc = options.cc
sq = options.sq
name = options.name
run = options.run

Cell_Concordance = pd.read_csv(cc,sep='\t')
Swap_Quant = pd.read_csv(sq,sep='\t')

Swap_Quant = Swap_Quant.set_index('cell')
Cell_Concordance = Cell_Concordance.set_index('GT 1')

Joined_Df = Swap_Quant.join(Cell_Concordance,how='inner')
Joined_Df['pool id']= name

try:
    cell_ambientness=pd.read_csv(f'/lustre/scratch123/hgi/projects/cardinal_analysis/qc/{run}/Donor_Quantification/{name}/ambientness_per_cell_{name}.tsv',sep='\t')
    cell_ambientness=cell_ambientness.set_index('barcode')
    Joined_Df = Joined_Df.join(cell_ambientness,how='inner')
except:
    _='No cell ambientless available'

Joined_Df = Joined_Df.reset_index()
Joined_Df.to_csv(f'{name}__joined_df_for_plots.tsv',sep='\t',index=False)



Joined_Df['total number of sites']=Joined_Df['Nr_Concordant']+Joined_Df['Nr_Discordant']

ax1 = sns.violinplot(data=Joined_Df, y="Percent_strict_discordant", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('becoming_different_donor.png')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('sites_becoming_different_donor.png')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="Total_reads", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('Total_reads_becoming_different_donor.png')
fig.clf()


Joined_Df2 = Joined_Df[Joined_Df["Nr times becoming different donor in subsampling"]!=0]
try:
    ax1 = sns.violinplot(data=Joined_Df2, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
    fig = ax1.get_figure()
    fig.savefig('sites_becoming_different_donor_no0.png')
    fig.clf()
except:
    _='There are no cells becoming different donor here.'

try:
    fig.clf()
    ax1 = sns.violinplot(data=Joined_Df2, y="Discordant_reads", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
    fig = ax1.get_figure()
    fig.savefig('Discordant_reads_becoming_different_donor_no0.png')
    fig.clf()
except:
    _='There are no cells becoming different donor here.'    

fig.clf()
ax1 = sns.violinplot(data=Joined_Df, y="Discordant_reads", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('Discordant_reads_becoming_different_donor.png')
fig.clf()

fig.clf()
ax1 = sns.violinplot(data=Joined_Df, y="Discordant_reads_by_n_sites", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('Discordant_reads_by_n_sites_becoming_different_donor.png')
fig.clf()


fig.clf()
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(13, 6))
sns.violinplot(data=Joined_Df, y="Nr_concordant_informative", x="Nr times becoming different donor in subsampling", cut=0, scale='width',ax=ax1)
sns.violinplot(data=Joined_Df, y="Nr_discordant_uninformative", x="Nr times becoming different donor in subsampling", cut=0, scale='width',ax=ax2)
fig.savefig('Nr_discordant_uninformative_becoming_different_donor.png')
fig.clf()

try:
    fig.clf()
    ax1 = sns.violinplot(data=Joined_Df2, y="Discordant_reads_by_n_sites", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
    fig = ax1.get_figure()
    fig.savefig('Discordant_reads_by_n_sites_becoming_different_donor_no0.png')
    fig.clf()
except:
    _='There are no cells becoming different donor here.'    

ax1 = sns.violinplot(data=Joined_Df, y="prob_max", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
# ax1 = sns.swarmplot(data=Joined_Df, y="prob_max", x="Nr times becoming different donor in subsampling",color= "white")
fig = ax1.get_figure()
fig.savefig('sites_becoming_different_donor_probs.png')
fig.clf()
ax1 = sns.violinplot(data=Joined_Df, y="Percent_strict_discordant", x="Nr times becoming Unassigned in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('becoming_unassigned_donor.png')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('sites_becoming_unassigned_donor.png')
fig.clf()

fig, ax1 = plt.subplots()
ax1 = sns.violinplot(data=Joined_Df, y="Percent_strict_discordant", x="Nr times becoming Doublet in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('becoming_doublet_donor.png')
ax1.hlines(y=0.2, xmin=0, xmax=20, linewidth=2, color='r')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('sites_becoming_doublet_donor.png')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="Total_reads", x="Nr times becoming different donor in subsampling", cut=0, scale='width')
fig = ax1.get_figure()
fig.savefig('Total_reads_becoming_different_donor.png')
fig.clf()


def scatter(fig, ax):
    
    sns.scatterplot(
        data=Joined_Df,
        x="Percent_strict_discordant",
        y="total number of sites",
        color="k",label=f"total nr cells assigned to donor={len(Joined_Df)}",
        ax=ax, alpha=0.5
    )

    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming Unassigned in subsampling']!=0]
    sns.scatterplot(
        data=Joined_Df_swap,
        x="Percent_strict_discordant",
        y="total number of sites",
        color="b",label=f"becoming Unassigned; total={len(Joined_Df_swap)}",
        ax=ax, alpha=0.5
    )
    # try:
    #     sns.kdeplot(
    #         data=Joined_Df_swap,
    #         x="Percent_strict_discordant",
    #         y="total number of sites",
    #         levels=2,
    #         color='b',
    #         fill=True,
    #         alpha=0.6,
    #         cut=2,
    #         ax=ax,
    #     )
    # except:
    #     _='only two entris, so cant do a density'

    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming Doublet in subsampling']!=0]
    sns.scatterplot(
        data=Joined_Df_swap,
        x="Percent_strict_discordant",
        y="total number of sites",
        color="y", label=f"becoming doublet; total={len(Joined_Df_swap)}",
        ax=ax, alpha=0.7
    )
    # try:
    #     sns.kdeplot(
    #         data=Joined_Df_swap,
    #         x="Percent_strict_discordant",
    #         y="total number of sites",
    #         levels=2,
    #         color='y',
    #         fill=True,
    #         alpha=0.6,
    #         cut=2,
    #         ax=ax,
    #     )
    # except:
    #     _='only two entris, so cant do a density'

    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']!=0]
    sns.scatterplot(
        data=Joined_Df_swap,
        x="Percent_strict_discordant",
        y="total number of sites",
        color="r", label=f"becoming different donor; total={len(Joined_Df_swap)}",
        ax=ax,
    )
    try:
        sns.kdeplot(
            data=Joined_Df_swap,
            x="Percent_strict_discordant",
            y="total number of sites",
            levels=3,
            color='r',
            fill=True,
            alpha=0.6,
            cut=2,
            ax=ax,
        )
    except:
        _='only two entris, so cant do a density'
        

    ax.legend()
    return fig

fig, ax = plt.subplots(figsize=(6, 6))
fig = scatter(fig, ax)
fig.savefig('sites_vs_concordance.png')
fig.clf()

try:
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(
        data=Joined_Df,
        y="Percent_strict_discordant",
        x="cell ambientness",
        color="k",label=f"total nr cells assigned to donor={len(Joined_Df)}",
        ax=ax, alpha=0.5
    )
    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']!=0]
    sns.scatterplot(
        data=Joined_Df_swap,
        y="Percent_strict_discordant",
        x="cell ambientness",
        color="r", label=f"becoming different donor; total={len(Joined_Df_swap)}",
        ax=ax,
    )

    fig.savefig('ambientness_vs_concordance.png')
    fig.clf()
except:
    _="Ambientness doesnt exist"

try:
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(
        data=Joined_Df,
        y="Discordant_reads",
        x="cell ambientness",
        color="k",label=f"total nr cells assigned to donor={len(Joined_Df)}",
        ax=ax, alpha=0.5
    )
    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']!=0]
    sns.scatterplot(
        data=Joined_Df_swap,
        y="Discordant_reads",
        x="cell ambientness",
        color="r", label=f"becoming different donor; total={len(Joined_Df_swap)}",
        ax=ax,
    )

    fig.savefig('ambientness_vs_read_concordance.png')
    fig.clf()
except:
    _="Ambientness doesnt exist"

try:
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(
        data=Joined_Df,
        y="Discordant_reads_by_n_sites",
        x="cell ambientness",
        color="k",label=f"total nr cells assigned to donor={len(Joined_Df)}",
        ax=ax, alpha=0.5
    )
    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']!=0]
    sns.scatterplot(
        data=Joined_Df_swap,
        y="Discordant_reads_by_n_sites",
        x="cell ambientness",
        color="r", label=f"becoming different donor; total={len(Joined_Df_swap)}",
        ax=ax,
    )

    fig.savefig('ambientness_vs_readbysites_concordance.png')
    fig.clf()
except:
    _="Ambientness doesnt exist"


import math
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Plot of plots, distinguish how many times the cell becomes a different donor
all_sub_times = set(Joined_Df['Nr times becoming different donor in subsampling'])
nr_plots=4
fig, axs = plt.subplots(math.ceil(len(all_sub_times)/2),nr_plots, figsize=(6*2, 6*math.ceil(len(all_sub_times)/2)),gridspec_kw={'width_ratios': [6, 1,6,1]})

st1=-1
st2=-2
all_sub_times=list(all_sub_times)
for i in range(len(all_sub_times)):
    print(i)
    i2=all_sub_times[i]
    if i % 2 == 0:
        # print(f"yes {i}")
        st1+=1
        # st2+=1
        st2=-2   
    st2+=2     
    try:
        ax1=axs[st1, st2]
        ax2=axs[st1, st2+1]
    except:
        ax1=axs[st2]
        ax2=axs[st2+1]

    # print(f"yes [{st1},{st2}]")
    # print(f"yes [{st1},{st2+1}]")
    if i==0:
        scatter(fig, ax1)
        try:
            Joined_Df_swap=Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']!=0]
            Joined_Df_swap['Nr times becoming different donor in subsampling']='at least once'
            sns.violinplot(data=Joined_Df_swap, y="total number of sites", x="Nr times becoming different donor in subsampling", scale='width',ax=ax2,color='r')
        except:
            _='no cells that swap donors'
        continue
    sns.scatterplot(
    data=Joined_Df,
    x="Percent_strict_discordant",
    y="total number of sites",
    color="k",
    ax=ax1, alpha=0.2
    )
    
    Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']==i2]
    
    try:
        sns.kdeplot(
            data=Joined_Df_swap,
            x="Percent_strict_discordant",
            y="total number of sites",
            levels=3,
            color='r',
            fill=True,
            alpha=0.6,
            cut=2,
            ax=ax1,
        )
    except:
        _='only two entris, so cant do a density'    
    sns.scatterplot(
        data=Joined_Df_swap,
        x="Percent_strict_discordant",
        y="total number of sites",
        color="r",
        ax=ax1, alpha=0.7,label=f"becoming different donor; n={i2}, total={len(Joined_Df_swap)}"
    )

    sns.violinplot(data=Joined_Df_swap, y="total number of sites", x="Nr times becoming different donor in subsampling", scale='width',ax=ax2,color='r')

        
ax.legend()
fig.tight_layout()
fig.savefig('subplot_sites_vs_concordance.png')
fig.clf()
print('Done')
