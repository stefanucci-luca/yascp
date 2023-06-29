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
    '-name', '--name',
    action='store',
    dest='name',
    required=True,
    help='name.'
)

options = parser.parse_args()

cc = options.cc
name = options.name

Cell_Concordance = pd.read_csv(cc,sep='\t')

Joined_Df = Cell_Concordance
Joined_Df['total number of sites']=Joined_Df['Nr_Concordant']+Joined_Df['Nr_Discordant']


ax1 = sns.violinplot(data=Joined_Df, y="Percent_strict_discordant", x="Nr times becoming different donor in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('becoming_different_donor.png')
fig.clf()

# Joined_Df = Joined_Df[Joined_Df["Nr times becoming different donor in subsampling"]!=0]
ax1 = sns.violinplot(data=Joined_Df, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('sites_becoming_different_donor.png')
fig.clf()

Joined_Df2 = Joined_Df[Joined_Df["Nr times becoming different donor in subsampling"]!=0]
ax1 = sns.violinplot(data=Joined_Df2, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('sites_becoming_different_donor_no0.png')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="Percent_strict_discordant", x="Nr times becoming Unassigned in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('becoming_unassigned_donor.png')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('sites_becoming_unassigned_donor.png')
fig.clf()

fig, ax1 = plt.subplots()
ax1 = sns.violinplot(data=Joined_Df, y="Percent_strict_discordant", x="Nr times becoming Doublet in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('becoming_doublet_donor.png')
ax1.hlines(y=0.2, xmin=0, xmax=20, linewidth=2, color='r')
fig.clf()

ax1 = sns.violinplot(data=Joined_Df, y="total number of sites", x="Nr times becoming different donor in subsampling", cut=0)
fig = ax1.get_figure()
fig.savefig('sites_becoming_doublet_donor.png')
fig.clf()

fig, ax = plt.subplots(figsize=(6, 6))
sns.scatterplot(
    data=Joined_Df,
    x="Percent_strict_discordant",
    y="total number of sites",
    color="k",
    ax=ax, alpha=0.5
)



Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming Unassigned in subsampling']!=0]
sns.scatterplot(
    data=Joined_Df_swap,
    x="Percent_strict_discordant",
    y="total number of sites",
    color="b",label="becoming Unassigned",
    ax=ax, alpha=0.5
)

sns.kdeplot(
    data=Joined_Df_swap,
    x="Percent_strict_discordant",
    y="total number of sites",
    levels=2,
    color='b',
    fill=True,
    alpha=0.6,
    cut=2,
    ax=ax,
)


Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming Doublet in subsampling']!=0]
sns.scatterplot(
    data=Joined_Df_swap,
    x="Percent_strict_discordant",
    y="total number of sites",
    color="y", label="becoming doublet",
    ax=ax, alpha=0.7
)

sns.kdeplot(
    data=Joined_Df_swap,
    x="Percent_strict_discordant",
    y="total number of sites",
    levels=2,
    color='y',
    fill=True,
    alpha=0.6,
    cut=2,
    ax=ax,
)

Joined_Df_swap = Joined_Df[Joined_Df['Nr times becoming different donor in subsampling']!=0]
sns.scatterplot(
    data=Joined_Df_swap,
    x="Percent_strict_discordant",
    y="total number of sites",
    color="r", label="becoming different donor",
    ax=ax,
)

sns.kdeplot(
    data=Joined_Df_swap,
    x="Percent_strict_discordant",
    y="total number of sites",
    levels=2,
    color='r',
    fill=True,
    alpha=0.6,
    cut=2,
    ax=ax,
)

ax.legend()
fig.savefig('sites_vs_concordance.png')

fig.clf()

print('Done')
