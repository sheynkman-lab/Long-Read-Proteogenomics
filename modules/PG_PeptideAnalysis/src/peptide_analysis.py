#!/usr/bin/env python3

""" 
This module prepares a table comparing mass spec MM peptide results using different databases

    Inputs:
    ------------------------------------------------------------------------------------------
    1. gene isoname file: map transcript name to gene name 
    2. Gencode peptides file: AllPeptides file from mass spec search using Gencode 
    3. Pacbio peptides file: Pacbio refined database fasta file 
    4. Pacbio six frame translation: file listing all possible peptides that can be detected per gene in Pacbio Database
    ------------------------------------------------------------------------------------------

    Output Tables:
    ------------------------------------------------------------------------------------------
    - table comparing pacbio coverage of Gencode peptide results from MM
    ------------------------------------------------------------------------------------------
"""
#%%
# Import Modules
import pandas as pd 
import re
import argparse
import os
from pathlib import Path 
from collections import defaultdict
from builtins import any  
from Bio import SeqIO

#%%
# Import Files
parser = argparse.ArgumentParser(description='Process peptide related input file locations')
parser.add_argument('--gene_to_isoname', '-gmap', action='store', dest='gene_isoname_file', help = 'Gene names to transcript names file location')
parser.add_argument('--gc_pep', '-gc', action='store', dest='gc_pep_file', help='Genecode AllPeptides file location')
parser.add_argument('--pb_pep', '-pb', action='store', dest='pb_ref_file', help='Pacbio AllPeptides file location')
parser.add_argument('--pb_6frm', '-sft', action='store', dest='pb_6frm_file', help='Pacbio Six Frame Translation file location')
parser.add_argument('--pb_gene', action='store', dest='pb_gene', help='PB to Gene file')
parser.add_argument("--cpat_all_orfs", action='store', dest='cpat_all_orfs')
parser.add_argument("--cpat_best_orf", action='store', dest='cpat_best_orf')
parser.add_argument("cpat_orf_protein_fasta", action='store', dest='cpat_protein_fasta')
parser.add_argument('-odir', '--output_directory', action='store', dest='odir', help = 'ouput directory')
results = parser.parse_args()

gene_isoname_file = results.gene_isoname_file
gc_pep_file = results.gc_pep_file
pb_refined_file = results.pb_ref_file 
pb_6frm_file = results.pb_6frm_file
pb_gene_file = results.pb_gene
cpat_best_orf_file = results.cpat_best_orf 
cpat_all_orfs_file = results.cpat_all_orfs 
cpat_protein_sequence_file = results.cpat_protein_fasta
#%%
# Input Filepaths 
gene_isoname_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/gene_isoname.tsv'
gc_pep_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/AllPeptides_Gencode.psmtsv'
pb_refined_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat_orf_aggregated.fasta'
pb_6frm_file ='/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.6frame.fasta'
pb_gene_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/pb_gene.tsv'
cpat_best_orf_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.ORF_prob.best.tsv'
cpat_protein_sequence_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.ORF_seqs.fa'
cpat_all_orfs_file = '/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/jurkat/jurkat.ORF_prob.tsv'
#%%

# loading gencode peptide data, initiate a dataframe
df = pd.read_table(gene_isoname_file, header=None)
isoname_gene = pd.Series(df[0].values, index=df[1]).to_dict()


# import gencode metamorpheus peptide data, filter to 1%FDR
g_cols = ['Base Sequence', 'Protein Accession', 'Decoy/Contaminant/Target', 'QValue']
g_data = pd.read_table(gc_pep_file, usecols = g_cols)
g_data.columns = ['pep_seq', 'acc', 'dct', 'qval']
g_tdata = g_data[(g_data['qval'] <= 0.01) & (g_data['dct']=='T')].reset_index(drop=True)
gc = g_tdata




# pb to gene map
df_pb_gene = pd.read_table(pb_gene_file)
pb_gene = pd.Series(df_pb_gene.gene.values, index=df_pb_gene.pb_acc).to_dict()

# replace each isoname with its gene name, explode distinct genes
def get_gene_name(row):
    isonames = re.split('\||\.(?=\D)', row['acc'])
    genes = set()
    for isoname in isonames:
        # TODO - fix the issue of unparsed isoanmes
        if isoname not in isoname_gene: continue # issues with lowercase parsing, Rob working on it 201122
        gene = isoname_gene[isoname]
        genes.add(gene)
    genes = list(genes)
    if len(genes) == 0:
        return 'no_match'
    return genes
gc['genes'] = gc.apply(get_gene_name, axis=1)



# TODO - debug, see TODO above
# print out isonames without a gene match
# found 282 peptides with no matched gene, Rob troubleshooting issue (with parsing of lowercase chars)
#gc[(gc['genes'] == 'no_match')]



# ~5K peptides duplicated in the allpeptides file, due to peptides with diff. mods identified
# gc[gc.duplicated(keep=False)]



gc = gc[['genes', 'pep_seq']]
gc = gc.explode('genes')
gc = gc.drop_duplicates()
gc.columns = ['gene', 'pep_seq']
gc_gene = gc.groupby('gene')['pep_seq'].apply(list).reset_index(name='gc_peps')
gc_gene['peps_in_gc'] = gc_gene['gc_peps'].apply(len)

# ~77K unique peptide-to-gene pairs
# 8018 unique genes

# presence of gc peptides in pb databse (generic function)
def get_pb_pep_coverage_stats(row, pb_dict):
    gene, peps, peps_in_gc = row['gene'], row['gc_peps'], row['peps_in_gc']
    if gene not in pb_dict:
        return 0, 0, 0
    else:
        num_peps_in_db = 0
        for pep in peps:
            if pep in pb_dict[gene]:
                num_peps_in_db += 1
        frac_peps_in_db = num_peps_in_db / peps_in_gc
        return 1, num_peps_in_db, frac_peps_in_db


## add in info for pb 6frm 
pb_6frm =  defaultdict()
for rec in SeqIO.parse(pb_6frm_file, 'fasta'):
    pb_6frm[rec.id] = str(rec.seq)
db = ('6frm', pb_6frm)
gc_gene[['in_{}'.format(db[0]),
         'peps_in_{}'.format(db[0]),
         'frac_peps_in_{}'.format(db[0])]] \
         = gc_gene.apply(lambda x: get_pb_pep_coverage_stats(x, db[1]), axis=1, result_type='expand')


## add in info for orf refined (ben) 
pb_refined = {}
for rec in SeqIO.parse(pb_refined_file, 'fasta'):
    gene = rec.description.split('=')[1]
    if gene not in pb_refined:
        pb_refined[gene] = ''
    pb_refined[gene] += '-' + str(rec.seq)

db = ('refined', pb_refined)
gc_gene[['in_{}'.format(db[0]),
         'peps_in_{}'.format(db[0]),
         'frac_peps_in_{}'.format(db[0])]] \
         = gc_gene.apply(lambda x: get_pb_pep_coverage_stats(x, db[1]), axis=1, result_type='expand')


#%%

## add in info for best cpat orfs
df_best_orfs = pd.read_table(cpat_best_orf_file)
#%%
df_best_orfs = df_best_orfs[['ID', 'Hexamer']]
df_best_orfs.columns = ['orf', 'coding_score']
df_best_orfs['pb_acc'] = df_best_orfs['orf'].str.split('_').str[0]
df_best_orfs['gene'] = df_best_orfs['pb_acc'].map(pb_gene)

#%%
# load in cpat prot seqs
pb_seq = defaultdict()
for rec in SeqIO.parse(cpat_protein_sequence_file, 'fasta'):
    pb_seq[rec.id] = str(rec.seq.translate())

# ...continuted cpat best orf
df_best_orfs['prot_seq'] = df_best_orfs['orf'].map(pb_seq)
df_best_orf_grp = df_best_orfs[['gene', 'prot_seq']].groupby('gene')['prot_seq'].apply(lambda x: '-'.join(x)).reset_index()
pb_best = pd.Series(df_best_orf_grp.prot_seq.values, index=df_best_orf_grp.gene).to_dict()
db = ('cpat_best', pb_best)
gc_gene[['in_{}'.format(db[0]),
         'peps_in_{}'.format(db[0]),
         'frac_peps_in_{}'.format(db[0])]] \
         = gc_gene.apply(lambda x: get_pb_pep_coverage_stats(x, db[1]), axis=1, result_type='expand')

## add in info for longest pb orf
cpat = pd.read_table(cpat_all_orfs_file)
cpat['pb_acc'] = cpat['ID'].str.split('_').str[0]
cpat = cpat.loc[cpat.groupby('pb_acc')['ORF'].idxmax()][['pb_acc', 'ID']]
cpat.columns = ['pb_acc', 'orf']
cpat['gene'] = cpat['pb_acc'].map(pb_gene)
cpat['prot_seq'] = cpat['orf'].map(pb_seq)
cpat = cpat[['gene', 'prot_seq']].groupby('gene')['prot_seq'].apply(lambda x: '-'.join(x)).reset_index()
pb_long = pd.Series(cpat.prot_seq.values, index=cpat.gene).to_dict()
db = ('cpat_long', pb_long)
gc_gene[['in_{}'.format(db[0]),
         'peps_in_{}'.format(db[0]),
         'frac_peps_in_{}'.format(db[0])]] \
         = gc_gene.apply(lambda x: get_pb_pep_coverage_stats(x, db[1]), axis=1, result_type='expand')

# If output directory DNE, make it
# odir = '../../results/PG_PeptideAnalysis'
odir = results.odir
if not os.path.exists(odir):
    os.mkdir(odir)

#%%
## write out file
gc_gene.to_csv(os.path.join(odir, 'gc_pb_overlap_peptides.tsv'), sep='\t', index=None)
# %%
