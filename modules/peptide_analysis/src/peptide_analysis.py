#!/usr/bin/env python3

""" 
This module prepares a table comparing mass spec MM peptide results from gencode 
against the fasta sequences of various orf calling methods

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
parser.add_argument('--gencode_peptides', '-gc', action='store', dest='gc_pep_file', help='Genecode AllPeptides file location')
parser.add_argument('--gene_to_isoname', '-gmap', action='store', dest='gene_isoname_file', help = 'Gene names to transcript names file location')
parser.add_argument('--pb_refined_fasta', action='store', dest='pb_refined_fasta', help='Pacbio refined database fasta file location')
parser.add_argument('--pb_filtered_fasta', action='store',dest='pb_filtered_fasta')
parser.add_argument('--pb_hybrid_fasta', action='store', dest='pb_hybrid_fasta')
parser.add_argument('--pb_gene', action='store', dest='pb_gene', help='PB to Gene file')
parser.add_argument('-odir', '--output_directory', action='store', dest='odir', help = 'ouput directory')
results = parser.parse_args()

gene_isoname_file = results.gene_isoname_file
gc_pep_file = results.gc_pep_file
pb_refined_file = results.pb_refined_fasta
pb_filtered_file = results.pb_filtered_fasta
pb_hybrid_file = results.pb_hybrid_fasta
pb_gene_file = results.pb_gene

#%%
# Input Filepaths
# results_dir = '/Users/bj8th/Documents/Sheynkman-Lab/Data/21-05-19_Jurkat' 
# name='jurkat'
# gene_isoname_file = f'{results_dir}/reference_tables/gene_isoname.tsv'
# gc_pep_file = f'{results_dir}/metamorpheus/gencode/search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv'
# pb_refined_file = f'{results_dir}/protein_gene_rename/{name}.protein_refined.fasta'
# pb_filtered_file = f'{results_dir}/protein_filter/{name}.filtered_protein.fasta'
# pb_hybrid_file = f'{results_dir}/hybrid_protein_database/{name}_hybrid.fasta'
# pb_gene_file = f'{results_dir}/transcriptome_summary/pb_gene.tsv'
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
gc_coverage = gc.groupby('gene')['pep_seq'].apply(list).reset_index(name='gc_peps')
gc_coverage['peps_in_gc'] = gc_coverage['gc_peps'].apply(len)

#%%
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

def add_pb_coverage_stats(fasta_file, gc_coverage, name):
    pb_seq = {}
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        gene = rec.description.split('=')[1]
        if gene not in pb_seq:
            pb_seq[gene] = ''
        pb_seq[gene] += '-' + str(rec.seq)
    for key, value in pb_seq.items():
        pb_seq[key] = value+'-'
    

    gc_coverage[[f'in_{name}',
         f'peps_in_{name}',
         f'frac_peps_in_{name}']] \
         = gc_coverage.apply(lambda x: get_pb_pep_coverage_stats(x, pb_seq), axis=1, result_type='expand')

add_pb_coverage_stats(pb_refined_file, gc_coverage, 'refined')
add_pb_coverage_stats(pb_filtered_file, gc_coverage, 'filtered')
add_pb_coverage_stats(pb_hybrid_file, gc_coverage, 'hybrid')
#%%



# If output directory DNE, make it
# odir = '../../results/PG_PeptideAnalysis'
odir = results.odir
if not os.path.exists(odir):
    os.mkdir(odir)

#%%
## write out file
gc_coverage.to_csv(os.path.join(odir, 'gc_pb_overlap_peptides.tsv'), sep='\t', index=None)
# %%
