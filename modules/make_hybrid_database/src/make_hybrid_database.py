#!/usr/bin/env python3
#%%
import pandas as pd 
from Bio import SeqIO
from collections import defaultdict
import argparse
#%%
parser = argparse.ArgumentParser()
parser.add_argument('--protein_classification',action='store',dest='pclass')
parser.add_argument('--gene_lens',action='store',dest='gene_lens')
parser.add_argument('--pb_fasta',action='store',dest='pb_fasta')
parser.add_argument('--gc_fasta',action='store',dest='gc_fasta')
parser.add_argument('--refined_info',action='store',dest='refined_info')
parser.add_argument('--pb_cds_gtf',action='store',dest='pb_cds_gtf')
parser.add_argument('--name',action='store',dest='name')
parser.add_argument('--lower_kb',action='store',dest='lower_kb',type=float,default=1)
parser.add_argument('--upper_kb',action='store',dest='upper_kb',type=float,default=4)
parser.add_argument('--lower_cpm',action='store',dest='lower_cpm',type=float,default=3)

args = parser.parse_args()

#%%
pclass = pd.read_table(args.pclass)

gene_lens = pd.read_table(args.gene_lens)

all_genes = set(pclass['pr_gene'])
# %%
lower_cpm = args.lower_cpm
lower_size = 1000*args.lower_kb
upper_size = 1000*args.upper_kb
#%%
cpm_genes = set(
    pclass
    .groupby('pr_gene')['CPM'].sum()
    .reset_index()
    .query(f'CPM>={lower_cpm}')['pr_gene']
)
    
#%%

size_genes = set(
    gene_lens
        .query(f'avg_len>={lower_size} and avg_len<={upper_size}')['gene'] 
)
goldilocks_genes = cpm_genes.intersection(size_genes)
#%%
pclass = pclass[pclass['pr_gene'].isin(goldilocks_genes)]
pclass_accs = set(pclass['pb'])
# %%

# %%
pb_fasta = defaultdict(list)
for record in SeqIO.parse(args.pb_fasta,'fasta'):
    gene = record.description.split('fullname GN=')[1].strip()
    if gene in goldilocks_genes:
        acc = record.id.split('|')[1].strip()
        if acc in pclass_accs:
            pb_fasta[gene].append(record)
#%%
gc_genes_all = set()
gc_fasta = defaultdict(list)
for record in SeqIO.parse(args.gc_fasta, 'fasta'):
    gene_split = record.description.split('GN=')
    gene = gene_split[1].strip()
    gc_genes_all.add(gene)
    if gene not in goldilocks_genes:
        gc_fasta[gene].append(record)
# %%
aggregated_fasta = []
for gene in gc_genes_all:
    if gene in gc_fasta:
        aggregated_fasta = aggregated_fasta + gc_fasta[gene]
    if gene in pb_fasta:
        aggregated_fasta = aggregated_fasta + pb_fasta[gene]

#%%
SeqIO.write(aggregated_fasta, f'{args.name}_hybrid.fasta', 'fasta')
# %%
# write aggreated refined info table - pb filtered 
refined = pd.read_table(args.refined_info)
refined = refined[refined['base_acc'].isin(pclass_accs)]
refined.to_csv(f'{args.name}_refined_high_confidence.tsv',sep='\t',index=False)

# write cds high confidence table
with open(args.pb_cds_gtf) as ifile, open(f'{args.name}_cds_high_confidence.gtf', 'w') as ofile:
    for line in ifile:
        pb_acc= line.split("|")[1]
        if pb_acc in pclass_accs:
            ofile.write(line)
# %%
# write high confidence genes to table
with open(f'{args.name}_high_confidence_genes.tsv', 'w') as ofile:
    ofile.write("GENE\n")
    for gene in goldilocks_genes:
        ofile.write(f'{gene}\n')