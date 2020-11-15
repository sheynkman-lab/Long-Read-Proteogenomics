# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
## Import Modules ##
#import weget # TODO: Add paths to MM outputs from wget
import numpy as np
import pandas as pd
from pathlib import Path
import csv
from collections import defaultdict

## Import Custom Modules ##
from m_gen_maps import GenMap

## Input Files ##
trans_to_gene_file = './trans_to_gene.tsv'
gtf_file = '../jurkat_analysis/a_gencode_gene_models/gencode.v35.annotation.gtf' # only need this if trans_to_gene file does not exist
pbacc_to_gene_file = './uniprot_acc_to_gencode_gene.tsv'

## MM Inputs ##
pep = './../Map_PacBio_Gene_Space/New_ForGloria/GENCODE_New_28Tryp/Task2-SearchTask/AllPeptides.psmtsv'
pb_pep = './../Map_PacBio_Gene_Space/ForGloria/28FractionBUTrypsin/PacBio_new/Task2-SearchTask/AllPeptides.psmtsv'
six_fa = './../Map_PacBio_Gene_Space/b_pacbio_6frm_database_gene_grouped.fasta/b_pacbio_6frm_database_gene_grouped.fasta'

## TODO's ##
# Fix uniprot gene mapping issues 
# Add path in GenMAp file to control output location 
# Put all of the analysis portion (including dictionaries?) into separate module 
# Start adding in code for peptide comparison


# %%
## Prepare Dictionaries for Genecode MM Outputs ##

# If Trans to Gene file does not exist, make it 
if Path(trans_to_gene_file).is_file()==False:
    # TODO: Add path in GenMap file to control output location
    GenMap(gtf_file, 'trans_to_gene')

# Make dictionary of transcript_name -> gene 
trans_to_gene = pd.read_csv(trans_to_gene_file, sep='\t')
trans_to_gene.columns = ['A', 'B']
gdict = pd.Series(trans_to_gene.A.values, index = trans_to_gene.B).to_dict()


# %%
# Import Peptides Dataset
cols = ['Base Sequence', 'Protein Accession', 'Decoy/Contaminant/Target', 'QValue']
data = pd.read_csv(pep, delimiter=r"\t", usecols = cols)
data.columns = ['seq', 'pb_acc', 'dct', 'qval']

# Filter Data
fdata = data[data['qval'] <= 0.01]
tdata = fdata[fdata['dct'] == 'T']

# Extract pb_acc col and split each row 
cut = tdata[['pb_acc']]
split = cut.pb_acc.str.split('\||\.', expand=True)

# Replace pb_acc (transcript_name) -> gene_name and remove duplicates
# TODO: There are some transcripts that don't map to genes. Ignoring those for now
gen = split.apply(lambda x: x.map(gdict, na_action='ignore')) 
tdata['gene'] = gen.stack().groupby(level=0).apply(lambda x: x.unique().tolist())


# Gene based results
p = tdata.gene.apply(pd.Series)
p.insert(0,'seq', tdata.seq.values)

# sort the sequences by gene
sort = p.melt(id_vars=['seq'], value_name="gene").dropna().reset_index(drop=True).drop('variable',1)
gen_pro = sort.groupby('gene')


# %%
# Import Peptides Dataset
cols = ['Base Sequence', 'Gene Name', 'Decoy/Contaminant/Target', 'QValue']
data = pd.read_csv(pb_pep, delimiter=r"\t", usecols = cols)
data.columns = ['seq', 'gene', 'dct', 'qval']

# Remove 'primary:' from gene column 
data['gene'] = data['gene'].str.replace('primary:', '')

# Filter Data
fdata = data[data['qval'] <= 0.01]
tdata = fdata[fdata['dct'] == 'T']

# Extract pb_acc col and split each row 
cut = tdata[['gene']]
split = cut.gene.str.split('\||\.', expand=True)

# Remove duplicate genes and sort by gene
tdata['gene'] = split.stack().groupby(level=0).apply(lambda x: x.unique().tolist())

# Gene based results
p = tdata.gene.apply(pd.Series)
p.insert(0,'seq', tdata.seq.values)

# sort the sequences by gene
sort = p.melt(id_vars=['seq'], value_name="gene").dropna().reset_index(drop=True).drop('variable',1)
pb_pro = sort.groupby('gene')


# %%
# Make dictionary for six frame translation of pacbio data
gene_seqs = defaultdict(lambda: set())

with open(six_fa) as f:
    for line in f:
        if line.startswith('>'):
            gene = line.split('>')[1]
            seq = next(f,'').split('-')
            seq = [i.strip() for i in seq]
            gene_seqs[gene].update(seq)


# %%
#gen_pro.apply(lambda x: x.equals(pb_pro.get_group(x.name)) if x.name in pb_pro.groups else False
for group in pb_pro.groups:
    if group in gen_pro.groups:
        A = pb_pro.get_group(group)
        B = gen_pro.get_group(group)
        C = pd.merge(A,B,how='inner', on='seq')
        D = len(C)
        E = len(A)


# %%
groups = []
gengroups_notin_pb = []
pb_match = []
gen_pos = []
score = []
for group in gen_pro.groups:
    if group in pb_pro.groups:
        pb_seq = pb_pro.get_group(group).drop_duplicates()
        gen_seq = gen_pro.get_group(group).drop_duplicates()
        merge_seq= pd.merge(pb_seq, gen_seq,how='inner', on='seq')
        gen_pep = len(gen_seq)
        pb_pep = len(merge_seq)
        pb_overlap = pb_pep/gen_pep
        # Save Values
        groups.append(group)
        pb_match.append(pb_pep)
        gen_pos.append(gen_pep)
        score.append(pb_overlap)
    if group not in pb_pro.groups:
        gen_seq = gen_pro.get_group(group).drop_duplicates()
        groups.append(group)
        gengroups_notin_pb.append(group)
        pb_match.append(0)
        gen_pos.append(len(gen_seq))
        score.append(0)


# %%
data = {'gene': groups, 'pb_pepmatch' : pb_match, 'gen_pep': gen_pos, 'pb_score' :score}
df = pd.DataFrame(data)
df.to_csv('./peptide_analysis.tsv', sep='\t', index=None)


# %%
import matplotlib.pyplot as plt

scores = df[['pb_score']]
scores.hist(bins=100)
#plt.hist(scores, bins=np.arange(min(scores), max(scores) + 0.01, 0.01))


# %%
len(gengroups_notin_pb)


# %%



