## Import Modules ##
# Python Modules #
import pandas as pd 
import csv
import argparse
from pathlib import Path 
from collections import defaultdict
from builtins import any  

# Custom Modules #
from gen_maps import GenMap 

# TODO: #
# Ask Gloria to put the Trans_to_Gene file on Zenodo #
# Figure out how to make optional command line arguments


# %%
## Import Files ##
""" 
Input Files:
-----------
- transcript to gene file: map transcript name to gene name
OR
- gtf file: contains jurkat genomic information 
- pbacc_to_gene file: map uniprot gene names to gencode gene names

Input Files from MM Outputs:
----------------------------
- Gencode Peptides File: AllPeptides file from genecode search 
- Pacbio Peptides File: AllPeptides file from pacbio search
- Pacbio Six Frame Translation: File listing all peptides that pacbio should be able to find 

Output Files:
-------------
- Table comparing pacbio coverage of gencode results
- TBD
"""
#parser = argparse.ArgumentParser(description='Process peptide related input file locations')
#parser.add_argument('--transcript_to_gene_file', '-tg', action='store', dest='trans_to_gene_file', help = 'File mapping transcript names to gene names location')
#parser.add_argument('--gtf_file', '-g', action='store', dest='gtf_file',help='GTF file location')
#parser.add_argument('--pbacc_to_gene', '-pbg', action='store', dest='pbacc_to_gene_file')
#parser.add_argument('--gen_pep', '-gp', action='store', dest='gen_pep_file', help='Genecode AllPeptides file location')
#parser.add_argument('--pb_pep', '-pp', action='store', dest='pb_pep_file', help='Pacbio AllPeptides file location')
#parser.add_argument('--pbsix_fram', '-sf', action='store', dest='six_fr_file', help='Pacbio Six Frame Translation file location')

## Input Files ##
trans_to_gene_file = '../../data/trans_to_gene.tsv'
gtf_file = '../../data/gencode.v35.annotation.gtf' # only need this if trans_to_gene file does not exist
pbacc_to_gene_file = '../../data/uniprot_acc_to_gencode_gene.tsv'

## MM Inputs ##
gen_pep_file = '../../data/AllPeptides_Gencode.psmtsv'
pb_pep_file = '../../data/AllPeptides_Pacbio.psmtsv'
six_fa = '../../data/pacbio_6frm_database_gene_grouped.fasta'


# %%
## Prepare Dictionaries for Genecode MM Outputs ##

# If trans to gene file does not exist, make it #
#TODO: Ask Gloria to put the trans_to_gene file in Zenodo 
if Path(trans_to_gene_file).is_file()==False:
    GenMap(gtf_file, 'trans_to_gene')

# Import File #
trans_to_gene = pd.read_csv(trans_to_gene_file, sep='\t')
trans_to_gene.columns = ['gene', 'trans']

# Map each transcript to a gene #
gencode_map = pd.Series(trans_to_gene.gene.values, index=trans_to_gene.trans).to_dict()


# %%
## Process Gencode Data ##

# Import AllPeptides File #
g_cols = ['Base Sequence', 'Protein Accession', 'Decoy/Contaminant/Target', 'QValue']
g_data = pd.read_csv(gen_pep_file, delimiter=r"\t", usecols = g_cols)
g_data.columns = ['seq', 'acc', 'dct', 'qval']

"After importing data we have 365553 peptides"

# Filter By Qvalue and Keep Targets #
g_fdata = g_data[g_data['qval'] <= 0.01]
g_tdata = g_fdata[g_fdata['dct']=='T'].reset_index(drop=True)

"After filtering we now have 76181 peptides"

"""
Now we have situations where each peptide is mapping to multiple accessions like this:
           GLATFCLDKDALRDEYDDLSDLNAVQMESVR        PGRMC2_201.PGRMC2_207|PGRMC2_208

So, I'm going to remove and split the accession column so it is like this instead:
                         PGRMC2_201 PGRMC2_207 PGRMC2_208
"""

# Extract acc col and split each row using . and | as delimiters #
g_cut = g_tdata[['acc']]
g_split = g_cut.acc.str.split('\||\.(?=\D)', expand=True)

# Replace each transcript acc with its gene name #
g_mapname = g_split.apply(lambda x: x.map(gencode_map, na_action='ignore'))

# Find locations where there was no gene_name match and fill them with "No_match" #
g_mapname.iloc[:,[0]] = g_mapname.iloc[:,[0]].fillna('No_match')

# Compile all of the gene names for each row and remove duplicates #
g_compile = g_mapname.stack().groupby(level=0).apply(lambda x: x.unique().tolist()).apply(pd.Series)

"We still have 76181 peptides"

# Add seqs, sort each peptide sequence by gene and then group by gene name #
g_compile.insert(0, 'seq', g_tdata.seq.values)
g_sort = g_compile.melt(id_vars=['seq'], value_name="gene").dropna().reset_index(drop=True).drop('variable',1)

"We now have 80389 rows, so there is less than 4208 cases where one peptide maps to more than 1 gene "

g_final = g_sort.groupby('gene')


# %%
## Process Pacbio Data ##

# Import Peptides Dataset #
cols = ['Base Sequence', 'Gene Name', 'Decoy/Contaminant/Target', 'QValue']
pb_data = pd.read_csv(pb_pep_file, delimiter=r"\t", usecols = cols)
pb_data.columns = ['seq', 'gene', 'dct', 'qval']

"After importing data we have 150806 peptides"

# Remove 'primary:' from gene column 
pb_data['gene'] = pb_data['gene'].str.replace('primary:', '')

# Filter By Qvalue and Keep Targets #
pb_fdata = pb_data[pb_data['qval'] <= 0.01]
pb_tdata = pb_fdata[pb_fdata['dct'] == 'T'].reset_index(drop=True)

"After filtering data we have 34660 peptides"

# Extract acc col and split each row using . and | as delimiters #
pb_cut = pb_tdata[['gene']]
pb_split = pb_cut.gene.str.split('\||\.(?=\D)', expand=True)

# Compile all of the gene names for each row and remove duplicates #
pb_compile = pb_split.stack().groupby(level=0).apply(lambda x: x.unique().tolist()).apply(pd.Series)

"We still have 34660 peptides"

# Add seqs, sort each peptide sequence by gene and then group by gene name #
pb_compile.insert(0, 'seq', pb_tdata.seq.values)
pb_sort = pb_compile.melt(id_vars=['seq'], value_name="gene").dropna().reset_index(drop=True).drop('variable',1)

"After sorting, we have 35271 rows, so there are less than 611 cases where one peptide maps to more than 1 gene"

pb_final = pb_sort.groupby('gene')


# %%
## Group Six Frame Translation Data ##

# Initialize Lists #
A_gene =[]
B_seq = []

# Parse file to get a list of genes and sequences #
with open(six_fa) as f:
    for line in f:
        if line.startswith('>'):
            gene = line.split('>')[1].strip()
            seq = next(f,'').replace(' ', '-').split('-')
            seq = [i.strip() for i in seq]
            for i in seq:
                A_gene.append(gene)
                B_seq.append(i)

# Convert lists to dataframe #
six_data = {'gene': A_gene, 'seq': B_seq}
PBsix = pd.DataFrame(six_data)
PBsix.columns = ['gene', 'seq']
PBsix = PBsix.sort_values(by=['gene', 'seq']).reset_index(drop=True)

# Group by Gene #
PBsix_final = PBsix.groupby('gene')


# %%
## Compare Peptides in Pacbio and Pacbio Sixframe Translation to Gencode ##

# Initialize Lists #
# For Pacbio Strict 
genes = []
pb_match = []
gen_num = []
score = []
pb_only = []
in_pb = []

# For Six Frame Translation 
s_genes = []
s_match = []
s_gen_num = []
s_score = []
six_only = []
in_six = []

# Main #
# For a gene in gencode
for group in g_final.groups:
    
    # If it is also in pacbio:
    if group in pb_final.groups:
        # Get peptides found in pacbio and genecode for a gene 
        pb_seq = pb_final.get_group(group).drop_duplicates()
        gen_seq = g_final.get_group(group).drop_duplicates()
        
        # Find overlapping genes
        union = pd.merge(pb_seq, gen_seq, how='inner', on='seq')
        overlap = len(union)/len(gen_seq)

        # Save Values to List
        genes.append(group)
        pb_match.append(len(union))
        gen_num.append(len(gen_seq))
        score.append(overlap)
        pb_only.append(overlap)
        in_pb.append('yes')

    if group not in pb_final.groups:
        # Get number of gencode peptides
        gen_seq = g_final.get_group(group).drop_duplicates()

        # Save Values to List
        genes.append(group)
        pb_match.append(0)
        gen_num.append(len(gen_seq))
        score.append(0)
        in_pb.append('no')
    
    if group in PBsix_final.groups:
        # Get list of peptides for gene
        six_seq = PBsix_final.get_group(group).drop_duplicates().seq.values.tolist()
        gen_seq = g_final.get_group(group).drop_duplicates()

        # Find Number of Overlapping Peptides
        six = 0 
        for i in gen_seq.seq:
            if any(i in x for x in six_seq):
                six = six + 1
        overlap = six/len(gen_seq)

        # Save Values to list 
        s_genes.append(group)
        s_match.append(six)
        s_gen_num.append(len(gen_seq))
        s_score.append(overlap)
        six_only.append(overlap)
        in_six.append('yes')
    
    if group not in PBsix_final.groups:
        # Get number of gencode peptides
        gen_seq = g_final.get_group(group).drop_duplicates()

        # Save Values to List
        s_genes.append(group)
        s_match.append(0)
        s_gen_num.append(len(gen_seq))
        s_score.append(0)
        in_six.append('no') 
        


# %%
# Compile Lists in DataFrame # 
data_pb = {'gene': genes, 'GC_pep': gen_num, 'in_pb': in_pb, 'PBorf_pep': pb_match, 'pb_overlap': score}
df1 = pd.DataFrame(data_pb)

data_sf = {'gene': s_genes, 'test_GC': s_gen_num, 'in_PB6frm': in_six, 'PB6frm_pep': s_match, 'PB6frm_overlap': s_score}
df2 = pd.DataFrame(data_sf)

data = pd.merge(df1, df2, how='outer', on='gene').drop(columns=['test_GC'])

#gen_tab = pd.read_csv('gene_based_info.tsv', sep='\t')
#tab = pd.merge(gen_tab, data, how='inner', on='gene')
#tab.to_csv('gencode_pb_comparison_tab.tsv', sep='\t', index=None)


# %%
# Plot Stuff #
import matplotlib.pyplot as plt

f = plt.figure(figsize=(20,15))
ax = f.add_subplot(2, 2, 1)
ax2 = f.add_subplot(2, 2, 2)
ax3 = f.add_subplot(2, 2, 3)
ax4 = f.add_subplot(2, 2, 4)
ax.hist(score, bins=100)
ax.set_title('PBorf_overlap')
ax2.hist(s_score, bins=100)
ax2.set_title('PB6frm_overlap')
ax3.hist(pb_only, bins=100)
ax3.set_title('PBorf_overlap for Gencode genes in Pacbio')
ax4.hist(six_only, bins=100)
ax4.set_title('PB6frm_overlap for Gencode genes in Pacbio6frm')
f.savefig('Pacbio_Coverage_of_Gencode_Genes.pdf')


