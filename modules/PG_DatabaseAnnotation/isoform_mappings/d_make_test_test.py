# make a subset of the isoform mapping file, only the toy genes


import pandas as pd

df = pd.read_table('c_uniprot_gencode_pacbio_isoform_mappings.tsv')


df.head()


toy_genes = [x.strip() for x in open('GenesForTestSet_Fraction16.txt').readlines()]
toy_genes


df.columns

df = df[['uniprot_acc', 'pacbio_acc', 'gencode_acc_y', 'gencode_acc_x', 'uniprot_gene_gc_matched', 'pacbio_gene', 'gencode_gene_x', 'gencode_gene_y', 'gene_match_cat_v2', 'gene_match_v2', 'gene_mismatches_v2']]
df.to_csv('d_isoform_and_gene_table_toy_genes.tsv', sep='\t', index=None)
