#%%
import pandas as pd
import numpy as np
from collections import defaultdict
import copy
import argparse
import gtfparse
import logging


sample_gtf = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.collapsed.gff'
agg_orfs = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat_orf_aggregated.tsv'
refined_orfs = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat_best_orf.tsv'
pb_gene = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/pb_gene.tsv'


representative_accessions = pd.read_table(agg_orfs)['base_acc'].to_list()
# import gtf, only exon info.
# only move forward with representative pb isoform (for same-protein groups)
gtf = gtfparse.read_gtf(sample_gtf)

gtf = gtf[['seqname', 'feature', 'start', 'end', 'strand', 'transcript_id']]
gtf = gtf[gtf['feature'] == 'exon']
gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
# only move forward with "base accession" (representative pb)

gtf = gtf[gtf.acc.isin(representative_accessions)]

# pb coords into dict
pbs = defaultdict(lambda: ['chr', 'strand', [], [], [],[]]) # pb -> [chr, strand, [start, end], [block lengths],[cum. block lengths], [prior cumulative block lengths]]
# PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
for i, row in gtf.iterrows():
    chr, feat, start, end, strand, acc = row
    pbs[acc][0] = chr
    pbs[acc][1] = strand
    pbs[acc][2].append([int(start), int(end)])
# sort all coords, calc blocks
for acc, infos in pbs.items():
    strand = infos[1]
    if strand == '+':
        infos[2] = sorted(infos[2])
    elif strand == '-':
        infos[2] = sorted(infos[2], reverse=True)
    infos[3] = np.array([end-start+1 for [start, end] in infos[2]])
    infos[4] = np.cumsum(infos[3])
    infos[5] = infos[4] - infos[3]


# read in the ranges of orf on pb transcripts
ranges = pd.read_table(refined_orfs)[['pb_acc', 'orf_start', 'orf_end', 'CPM']]
ranges = ranges[ranges['pb_acc'].isin(representative_accessions)]

# read in pb to genename
pb_gene = pd.read_table(pb_gene)
pb_gene = pd.Series(pb_gene.gene.values, index=pb_gene.pb_acc).to_dict()

#%%
from Bio import SeqIO
collapsed_fasta_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.collapsed.fasta'
seqs = {}
for rec in SeqIO.parse(collapsed_fasta_file, "fasta"):
    pb_acc = rec.id.split('|')[0]
    seqs[pb_acc] = len(rec.seq)


# %%
with open("small_gtf.tsv", "w") as small, open("bigger_gtf.tsv","w") as big:
    small.write("acc\tgtf_size\tseq_size\n")
    big.write("acc\tgtf_size\tseq_size\n")
    for acc, infos in pbs.items():
        if str(acc).startswith("PB."):
            # gtf_size = infos[4][-1]
            gtf_size = np.sum(infos[3])
            seq_size = seqs[acc]

            if gtf_size < seq_size:
                small.write(f"{acc}\t{gtf_size}\t{seq_size}\n")
            if gtf_size > seq_size:
                big.write(f"{acc}\t{gtf_size}\t{seq_size}\n")
# %%
small = pd.read_table("small_gtf.tsv")
big = pd.read_table("bigger_gtf.tsv")
#%%
print(len(small))
print(len(big))

gtf_total = len(pbs)
combined = len(big) + len(small)

print(combined)
print(gtf_total)
print(combined/gtf_total)

print(len(small)/ gtf_total)
print(len(big)/ gtf_total)
# %%

# %%
small['delta'] = small['gtf_size'] - small['seq_size']
big['delta'] = big['gtf_size'] - big['seq_size']
both = pd.concat([small, big])
# %%
import matplotlib.pyplot as plt
both['delta'].hist(bins = 50, ranges=(-10,10))
# %%
both['rel_delta'] = both['delta'] / both['seq_size']

# %%
both['rel_delta'].hist(bins = 50)
# %%
both = pd.merge(both, ranges, how = 'inner', left_on = 'acc', right_on='pb_acc')
# %%
both['CPM-bin'] = pd.cut(both['CPM'], 10)

# %%
import seaborn as sns
sns.violinplot(data=both, x='CPM-bin', y = 'delta')


# %%
plt.scatter(x = both['CPM'], y=both['delta'], alpha = 0.05)
plt.xlim(0, 30)
plt.ylim(-400, 100)
# %%
sns.regplot(data=both, x= 'CPM', y = 'delta',scatter_kws={'alpha':0.05})
plt.xlim(0, 30)
plt.ylim(-400, 100)
# %%
sns.pairplot(both)
# %%
from Bio import SeqIO
collapsed_fasta_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.collapsed.fasta'
seqs_seq = {}
for rec in SeqIO.parse(collapsed_fasta_file, "fasta"):
    pb_acc = rec.id.split('|')[0]
    seqs_seq[pb_acc] = str(rec.seq)

# %%
