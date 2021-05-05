#!/usr/bin/env python3

# Find novel peptides that have been detected from a MetaMorpheus
# MS search against the PacBio database.


#%%

import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--pacbio_peptides',action='store',dest='peptides')
parser.add_argument('--gencode_fasta',action='store',dest='gencode_fasta')
parser.add_argument('--name',action='store',dest='name')
args = parser.parse_args()

# read in pacbio peptides
peps = pd.read_table(args.peptides)
peps = peps[(peps['QValue']<=0.01) & (peps['Decoy/Contaminant/Target']=='T')]
peps = peps[['Gene Name', 'Protein Accession', 'Base Sequence', 'Score', 'QValue']]
peps.columns = ['gene', 'acc', 'seq', 'score', 'qval']

def get_first_gene_name(gene_str):
    gene = gene_str.split(':')[1]
    if '|' in gene:
        gene = gene.split('|')[0]
    return gene

peps['gene'] = peps['gene'].apply(get_first_gene_name)

# all detected peptides (pacbio sample-specific)
peps_pb = peps[peps['acc'].str.startswith('PB')]
sample_peptides = peps_pb['seq'].to_list()

# import all the sequences from gencode into a big string
gc = [str(rec.seq) for rec in SeqIO.parse(args.gencode_fasta, 'fasta')]
gc_agg = ','.join(gc)

# find novel peptides
pacbio_peps = []
for pep in sample_peptides:
    # some peptides are indistinguishable (have I/L), take first one
    if '|' in pep:
        pep = pep.split('|')[0]
    # is the base peptide sequence in the ref database
    if pep not in gc_agg:
        pacbio_peps.append(pep)

# write out the pacbio accession and genename for each novel peptide
peps_novel = peps[peps['seq'].isin(pacbio_peps)]
peps_novel.to_csv(f'{args.name}.pacbio_novel_peptides.tsv', sep='\t', index=None)

