# derive a same-protein-sequence gencode file

# %%

import pandas as pd
from collections import defaultdict
from Bio import SeqIO


# make gencode clusters of same-protein-sequence entries
# evaluate protein sequence similarity within the same gene
gc = defaultdict(lambda: defaultdict(list)) # prot_seq -> gene -> [isonames]
for rec in SeqIO.parse('../../data/gencode.v35.pc_translations.fa', 'fasta'):
    gene = rec.id.split('|')[6]
    isoname = rec.id.split('|')[5]
    gc[rec.seq][gene].append(isoname)


# write a fasta, with same-protein sequences clustered
# for each group of same-protein-sequence isonames, choose the one 
with open('gencode.fasta', 'w') as ofile, open('gencode_isoname_clusters.tsv', 'w') as ofile2:
    ofile2.write('repr_isoname\tsame_prot_seq_isonames\n')
    for seq in gc:
        for gene, isonames in gc[seq].items():
            isonames = sorted(isonames)
            repr_isoname = isonames[0]
            ofile.write('>gc|{}| GN={}\n{}\n'.format(repr_isoname, gene, seq))
            ofile2.write(repr_isoname + '\t' + ','.join(isonames[1:]) + '\n')

