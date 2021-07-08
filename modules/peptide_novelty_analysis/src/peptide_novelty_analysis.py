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
parser.add_argument('--uniprot_fasta', action='store',dest='uniprot_fasta')
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
gencode_seq = [str(rec.seq) for rec in SeqIO.parse(args.gencode_fasta, 'fasta')]
gencode_all_sequences = ','.join(gencode_seq)

uniprot_seq = [str(rec.seq) for rec in SeqIO.parse(args.uniprot_fasta, 'fasta')]
uniprot_all_sequences = ','.join(uniprot_seq)
# find novel peptides
novel_peps = set()
novel_peps_to_gencode = set()
novel_peps_to_uniprot = set()

for pep in sample_peptides:
    # some peptides are indistinguishable (have I/L), take first one
    if '|' in pep:
        pep = pep.split('|')[0]
    # is the base peptide sequence in the ref database
    if pep not in gencode_all_sequences:
        novel_peps_to_gencode.add(pep)
    if pep not in uniprot_all_sequences:
        novel_peps_to_uniprot.add(pep)
novel_peps = novel_peps_to_gencode.intersection(novel_peps_to_uniprot)

# write out the pacbio accession and genename for each novel peptide
peps_novel = peps[peps['seq'].isin(novel_peps)]
peps_novel.to_csv(f'{args.name}.pacbio_novel_peptides.tsv', sep='\t', index=None)

peps_novel_to_gencode = peps[peps['seq'].isin(novel_peps_to_gencode)]
peps_novel_to_gencode.to_csv(f'{args.name}.pacbio_novel_peptides_to_gencode.tsv', sep='\t', index=None)


peps_novel_to_uniprot= peps[peps['seq'].isin(novel_peps_to_uniprot)]
peps_novel_to_uniprot.to_csv(f'{args.name}.pacbio_novel_peptides_to_uniprot.tsv', sep='\t', index=None)
