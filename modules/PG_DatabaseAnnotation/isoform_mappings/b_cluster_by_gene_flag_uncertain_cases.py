# cluster the isoforms to their mother genes
# determine which mapped isoforms may need to be removed


import pandas as pd


# derive genes for each accession mapping

df = pd.read_table('a_pacbio_gencode_uniprot_consolid_table.tsv')

# Step 1 - extract genes for uniprot and gencode
df['uniprot_gene_orig'] = df['uniprot_acc'].str.split('|').str[2].str.split('_HUMAN').str[0]
df['gencode_gene_x'] = df['gencode_acc_x'].str.split('|').str[6]
df['gencode_gene_y'] = df['gencode_acc_y'].str.split('|').str[6]


# Step 2 - get genes for pacbio

# Step 2.A derive pacbio to gene map
iso_ensg = pd.read_table('../../../a_SQANTI3_out/jurkat_classification.txt')[['isoform', 'associated_gene']]
iso_ensg.columns = ['pb_acc', 'ensg']
genes = pd.read_table('../../../a_gencode_gene_models/b_gencode_gene_maps/a_ensg_gene.tsv', header=None)
genes.columns= ['ensg', 'gene']
tmp = iso_ensg.merge(genes, on='ensg', how='left')
iso_gene = pd.Series(tmp.gene.values, index=tmp.pb_acc.values)


# Step 2.B add in pacbio genes
# if more than one gene corresponds, list them all

def get_associated_genes_as_string(pb_accs):
    genes = set()
    for acc in pb_accs:
        genes.add(iso_gene[acc])
    for gene in genes:
        if pd.isnull(gene):
            return ''
    else:
        genes = sorted(list(genes))
        gene_str = '|'.join(genes)
        return gene_str

def retrieve_all_associated_genes(pacbio_acc_str):
    if pd.isnull(pacbio_acc_str):
        return ''
    else:
        pb_accs = pacbio_acc_str.split('|')
        gene_str = get_associated_genes_as_string(pb_accs)
        return gene_str

df['pacbio_gene'] = df['pacbio_acc'].apply(retrieve_all_associated_genes)



# Step 3 - determine if genes are concordant and mark differing ones


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

df[['gene_match_cat', 'gene_match', 'gene_mismatches']] = df[['uniprot_gene', 'gencode_gene_x', 'gencode_gene_y', 'pacbio_gene']].apply(determine_gene_match_status, axis=1)

df.to_csv('b_isoform_mappings_w_genes.tsv', sep='\t', index=None)

# %%

# determine number of genes that have isoform mapping issues
all_match_genes = set()
pb_included_genes = set()
map_issue_genes = set()
map_issue_genes_only_pb = set()
for i, row in df.iterrows():
    if not pd.isnull(row.gene_mismatches):
        mismatched_genes = row.gene_mismatches.split('|')
        for g in mismatched_genes:
            map_issue_genes.add(g)
    if not pd.isnull(row.gene_mismatches) and not pd.isnull(row.pacbio_acc):
        mismatched_genes = row.gene_mismatches.split('|')
        for g in mismatched_genes:
            map_issue_genes_only_pb.add(g)
    if not pd.isnull(row.gene_match):
        all_match_genes.add(row.gene_match)
    if not pd.isnull(row.pacbio_acc) and not pd.isnull(row.gene_match):
        pb_included_genes.add(row.gene_match)

# %%

print(len(all_match_genes))
print(len(pb_included_genes))
print(len(map_issue_genes))
print(len(map_issue_genes_only_pb))
print(len(pb_included_genes - map_issue_genes))
print(len(pb_included_genes - map_issue_genes_only_pb))

# %%
with open('tmp.tsv', 'w') as ofile:
    ofile.write(str(map_issue_genes))
