# harmonize the uniprot/gencode genenames

import pandas as pd

df = pd.read_table('b_isoform_mappings_w_genes.tsv')


# uniprot to gencode gene mappings
un_gc = pd.read_table('../d_uniprot_to_gencode/c_ucsc_uniprot_track_aln_file/a_uniprot_acc_to_gencode_gene.tsv')[['uniprot_acc', 'gencode_gene']]
un_gc_dict = pd.Series(un_gc['gencode_gene'].values, index=un_gc['uniprot_acc']).to_dict()


def get_gencode_gene_to_which_uniprot_acc_maps(un_acc):
    if pd.isnull(un_acc):
        return ''
    un_acc_short = un_acc.split('|')[1]
    if un_acc_short in un_gc_dict:
        return un_gc_dict[un_acc_short]
    else:
        return ''

df['uniprot_gene_gc_matched'] = df['uniprot_acc'].apply(get_gencode_gene_to_which_uniprot_acc_maps)


def determine_gene_match_status(genes):
    filtered_genes = set()
    # do not include nan or blank genes
    for gene in genes:
        if pd.isnull(gene): continue
        if gene == '': continue
        filtered_genes.add(gene)
    if len(filtered_genes) == 1:
        return pd.Series(['match', list(filtered_genes)[0], ''])
    else:
        gene_str = '|'.join(sorted(list(filtered_genes)))
        return pd.Series(['mismatch', '', gene_str])

df[['gene_match_cat_v2', 'gene_match_v2', 'gene_mismatches_v2']] = df[['uniprot_gene_gc_matched', 'gencode_gene_x', 'gencode_gene_y', 'pacbio_gene']].apply(determine_gene_match_status, axis=1)

# %%

df.to_csv('c_uniprot_gencode_pacbio_isoform_mappings.tsv', sep='\t', index=None)
