import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json

################################################################################################
mirgene_gff_path = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\raw\hsa.gff'
mirbase_gff_path = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\raw\mirbase_hsa.gff'

figure_out = r'C:\Users\felix\PycharmProjects\human_mirna_evol\figures\mirgene_vs_mirbase_lengthdist.png'
################################################################################################


def read_mirgene(path):
    outdir = {'type': [], 'length': [], 'source': []}
    with open(path, 'r') as fh:
        for line in fh:
            if not line.startswith('#'):
                linedata = line.strip().split()
                type = linedata[2]
                length = int(linedata[4]) - int(linedata[3])
                outdir['type'].append(type)
                outdir['length'].append(length)
                outdir['source'].append('MirGeneDB')
    return outdir


def read_mirbase(path):
    outdir = {'type': [], 'length': [], 'source': []}
    with open(path, 'r') as fh:
        for line in fh:
            if not line.startswith('#'):
                linedata = line.strip().split()
                if linedata[2] == 'miRNA_primary_transcript':
                    type = 'pre_miRNA'
                elif linedata[2] == 'miRNA':
                    type = 'miRNA'
                length = int(linedata[4]) - int(linedata[3])
                outdir['type'].append(type)
                outdir['length'].append(length)
                outdir['source'].append('miRBase')
    return outdir


################################################################################################

mirgene = read_mirgene(mirgene_gff_path)
mirgene_df = pd.DataFrame.from_dict(mirgene)

mirbase = read_mirbase(mirbase_gff_path)
mirbase_df = pd.DataFrame.from_dict(mirbase)


df = pd.concat([mirgene_df, mirbase_df])

pre_df = df[df['type'] == 'pre_miRNA']
mirna_df = df[df['type'] == 'miRNA']

fig, axs = plt.subplots(1, 2)

l1 = sns.kdeplot(data=pre_df, x='length', hue='source', ax=axs[0], bw_adjust=2)
l2 = sns.kdeplot(data=mirna_df, x='length', hue='source', bw_adjust=2, ax=axs[1])
axs[0].set_title('pre-miRNAs')
axs[1].set_title('mature and star')
plt.savefig(figure_out)
plt.show()
print(df)

mirgenedb_pre_median = mirgene_df['length'][mirgene_df['type'] == 'pre_miRNA'].median()
print(mirgenedb_pre_median)





# sns.kdeplot(data=mirgene_df, x='length', hue='type', bw_adjust=2)
# plt.figure()
# sns.kdeplot(data=mirbase_df, x='length', hue='type', bw_adjust=2)
# plt.show()