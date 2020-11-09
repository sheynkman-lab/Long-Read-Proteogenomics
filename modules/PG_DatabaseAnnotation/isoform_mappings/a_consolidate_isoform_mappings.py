# consolidate all the mappings into a master table

import pandas as pd
import numpy as np
from Bio import SeqIO

pb_gc = pd.read_table('../cb_blast_results_copied_from_rivanna/pacbio_at_len_matches_to_gencode.tsv')
pb_un = pd.read_table('../cb_blast_results_copied_from_rivanna/pacbio_at_len_matches_to_uniprot.tsv')
un_gc = pd.read_table('../cb_blast_results_copied_from_rivanna/uniprot_at_len_matches_to_gencode.tsv')


print(len(pb_gc.gencode_acc))
print(len(set(pb_gc.gencode_acc)))

pg = pb_gc[['pacbio_acc', 'gencode_acc']]
pu = pb_un[['pacbio_acc', 'uniprot_acc']]
ug = un_gc[['uniprot_acc', 'gencode_acc']]

df = pd.merge(pg, pu, how='outer', on='pacbio_acc')
df2 = pd.merge(df, ug, how='outer', on='uniprot_acc')
df2 = df2.sort_values('pacbio_acc')
conditions = [(df2['gencode_acc_x'] == '') | (df2['gencode_acc_y'] == ''), df2['gencode_acc_x'] == df2['gencode_acc_y'], df2['gencode_acc_x'] != df2['gencode_acc_y']]
choices = ['-', '1', '0']
df2['gencode_match'] = np.select(conditions, choices)

# add in pb-only, gc-only, and un-only isoform accessions
pb_accs = pd.read_table('../../a_extract_and_cluster_pacbio_prot_seqs/a_pacbio_protein_clusters.tsv')['pb_accs'].to_list()
pb_df = pd.DataFrame(pb_accs, columns = ['pacbio_acc'])
un_accs = [rec.id for rec in SeqIO.parse('../../../a_protein_databases/uniprot_reviewed_canonical_and_isoform.fasta', 'fasta')]
un_df = pd.DataFrame(un_accs, columns = ['uniprot_acc'])
gc_accs = [rec.id for rec in SeqIO.parse('../b_gencode_fasta_same_prot_clustered.v35.fa', 'fasta')]
gc_df = pd.DataFrame(gc_accs, columns = ['gencode_acc'])
# series of joins to add pb/gc/un-only rows
df2 = pd.merge(pb_df, df2, on='pacbio_acc', how='outer')
df2 = pd.merge(gc_df, df2, left_on='gencode_acc', right_on='gencode_acc_x', how='outer')
df2 = pd.merge(un_df, df2, on='uniprot_acc', how='outer')


df2['pacbio_acc_dups'] = df2.duplicated(subset='pacbio_acc', keep=False)*1
df2['gencode_acc_x_dups'] = df2.duplicated(subset='gencode_acc_x', keep=False)*1
df2['uniprot_acc_dups'] = df2.duplicated(subset='uniprot_acc', keep=False)*1

df2.to_csv('a_pacbio_gencode_uniprot_consolid_table.tsv', sep='\t', index=None)
