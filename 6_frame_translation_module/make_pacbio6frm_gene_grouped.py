# make a six frame translation database for from the pacbio transcripts

from Bio import SeqIO
import pandas as pd
from collections import defaultdict

# inputs - modify
sqanti_iso_annot_fpath = '../SQANTI3_out/jurkat_classification.txt'
ensg_gene_table_fpath = '../ensg_gene.tsv'
pacbio_transcript_fasta_fpath = '../jurkat_corrected.fasta'

# get associated gene for each pb acc, for final write-out (below)
# pb_gene is pb_acc -> gene dictionary
df = pd.read_table(sqanti_iso_annot_fpath)[['isoform', 'associated_gene']]
df2 = pd.read_table(ensg_gene_table_fpath, header=None)
df2.columns = ['associated_gene', 'gene']
df3 = pd.merge(df, df2, on='associated_gene', how='left').fillna('NOVEL')[['isoform', 'gene']]
pb_gene = pd.Series(df3.gene.values, index=df3.isoform).to_dict()


# read in pacbio transcripts

def return_all_orfs(translated_seq):
    # return all contiguous polypeptide sequences 7 AA or higher
    seqs = translated_seq.split('*')
    filtered_seqs = []
    for seq in seqs:
        if len(seq) >= 7:
            filtered_seqs.append(str(seq))
    return filtered_seqs

gene_seqs = defaultdict(lambda: set()) # gene -> pacbio sequences as list

for rec in SeqIO.parse(pacbio_transcript_fasta_fpath, 'fasta'):
    gene = pb_gene[rec.id]
    F1 = rec.seq.translate()
    F2 = rec.seq[1:].translate()
    F3 = rec.seq[2:].translate()
    R1 = rec.seq.reverse_complement().translate()
    R2 = rec.seq.reverse_complement()[1:].translate()
    R3 = rec.seq.reverse_complement()[2:].translate()
    translations = [F1, F2, F3, R1, R2, R3]
    for tr in translations:
        orfs = set(return_all_orfs(tr))
        gene_seqs[gene].update(orfs)


# write out fasta file in which entries represent each gene and the
# pacbio-derived protein "space"
with open('pacbio_6frm_database_gene_grouped.fasta', 'w') as ofile:
    for gene, orfs in gene_seqs.items():
        ofile.write('>' + gene + '\n' + '-'.join(orfs) + '\n')
