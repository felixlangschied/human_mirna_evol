import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json

###################################################################################################################
snps = '/home/felixl/PycharmProjects/human_miRNA_evol/data/processed/SNPs/no_snps_per_region.json'
maf_path = '/home/felixl/PycharmProjects/human_miRNA_evol/data/processed/SNPs/mafs_per_region.json'
figout = '/home/felixl/PycharmProjects/human_miRNA_evol/figures/snp_density_mirnas.png'

###################################################################################################################

def load_snps(path):
    with open(path, 'r') as fh:
        tmpdic = json.load(fh)
    df = pd.DataFrame.from_dict(tmpdic, orient='index')
    df = df.drop(index='pre_mirna')
    df['ratio'] = df['numsnps'] / df['length'] * 1000
    return df


def categorize_mafs(maf_dict):
    common = 0
    uncommon = 0
    rare = 0
    ultra_rare = 0
    no_data = 0
    num_snps = len(maf_dict)

    for maf in maf_dict:
        try:
            maf = float(maf)
        except ValueError:
            try:
                maf = float(maf.split(';')[0])
            except:
                no_data += 1
                continue
        if maf >= 0.05:
            common += 1
        elif maf < 0.05 and maf >= 0.01:
            uncommon += 1
        elif maf < 0.01 and maf >= 0.0001:
            rare += 1
        elif maf < 0.0001:
            ultra_rare += 1
    return map(lambda x: x / num_snps, [common, uncommon, rare, ultra_rare, no_data])


def load_mafs(path):
    df_dict = {'count': [], 'rarity': [], 'region': []}
    with open(path, 'r') as fh:
        m_dict = json.load(fh)
        for region in m_dict:
            if region == 'pre_mirna':
                continue
            if region != 'mature':
                continue

            common, uncommon, rare, ultra_rare, no_data = categorize_mafs(m_dict[region])
            catlist = ['common', 'uncommon', 'rare', 'ultra_rare', 'no_data']

            df_dict['region'].append(region)
            df_dict['count'].append(common)
            df_dict['rarity'].append('common')

            df_dict['region'].append(region)
            df_dict['count'].append(uncommon)
            df_dict['rarity'].append('uncommon')

            df_dict['region'].append(region)
            df_dict['count'].append(rare)
            df_dict['rarity'].append('rare')

            df_dict['region'].append(region)
            df_dict['count'].append(ultra_rare)
            df_dict['rarity'].append('ultra_rare')

            df_dict['region'].append(region)
            df_dict['count'].append(no_data)
            df_dict['rarity'].append('no_data')
            break

    df = pd.DataFrame.from_dict(df_dict)
    return df

###########################################################################################################
### SNPs
###########################################################################################################
snpdf = load_snps(snps)


sns.set_style('whitegrid')
sns.set_theme('paper')
order = ['downstream', 'mature', 'star', 'hairpin', 'upstream', 'protein_coding', 'lncRNA']
ax = sns.barplot(data=snpdf, x=snpdf.index, y='ratio', order=order)
plt.ylabel('SNP density per 1 kb')

xticks = ['5p-Flank', 'Mature', 'Star', 'Hairpin', '3p-Flank', 'Protein coding', 'lncRNA']
ax.set_xticklabels(xticks, rotation=20, ha='right')
plt.tight_layout()
plt.savefig(figout, dpi=600)


######################################################################################################
### MAF
######################################################################################################

plt.figure()
mafdf = load_mafs(maf_path)
sns.barplot(data=mafdf, x='rarity', y='count', hue='region')
plt.show()