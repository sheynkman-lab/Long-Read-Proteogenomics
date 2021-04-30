#!/usr/bin/env python3

""" 
This module prepares reference tables for other modules 
 
  Inputs:
  --------------------------------------------------------------
   1. Gencode gtf file
   2. Gencode fasta file 
  --------------------------------------------------------------
  
   Output Tables:
  -------------------------------------------------------------
   1. ensg -> gene 
   2. isoname (transcript name) -> gene 
   3. ensp -> gene 
   4. isoform, gene, length table 
   5. gene -> min, max, average transcript length
   6. protein coding genes
 -------------------------------------------------------------- 

"""

# Import Modules 
import pandas as pd 
import argparse
import csv
import os
from collections import defaultdict
import gtfparse

# Define Functions
def GenMap(gtf_file, results):
    """
    This function prepares a series of tables that map gene, transcript and protein level information. 
    """

    # Initialize Lists and Dictionaries
    genes = {} # ENSG -> <genename>
    isos = {} # ENST -> <iso-name>
    ensps = defaultdict(set) # gene_name -> <set(ENSPs)>
    isonames = defaultdict(set) # transcript name -> gene_name

    # Parse GTF file
    for line in open(gtf_file):
        if line.startswith('#'):
            pass
        else:
            wds = line.split('\t')
            cat = wds[2]
            if cat in ['transcript']:
                ensg = line.split('gene_id "')[1].split('"')[0]
                gene = line.split('gene_name "')[1].split('"')[0]
                enst = line.split('transcript_id "')[1].split('"')[0]
                transcript_name = line.split('transcript_name "')[1].split('"')[0]
                genes[ensg] = gene
                isos[enst] = transcript_name
                isonames[gene].add(transcript_name)
                if 'transcript_type "protein_coding' in line:
                        gen = line.split('gene_name "')[1].split('"')[0]
                        ensp = line.split('protein_id "')[1].split('"')[0]
                        ensps[gen].add(ensp)
    
    # Save Tables in results/PG_ReferenceTables 
    with open(results.ensg_gene, 'w') as ofile:
        for ensg, gene in genes.items():
            ofile.write(ensg + '\t' + gene + '\n')
    
    with open(results.enst_isoname, 'w') as ofile:
        for enst, isoname in isos.items():
            ofile.write(enst + '\t' + isoname + '\n')
    
    with open(results.gene_ensp, 'w') as ofile:
        for gen, ensp_set in ensps.items():
            for ensp in ensp_set:
                ofile.write(gen + '\t' + ensp + '\n')
    
    with open(results.gene_isoname, 'w') as ofile:
        for gene, transcript_set in isonames.items():
            for transcript in transcript_set:
                ofile.write(gene + '\t' + transcript + '\n')
    
    print("The ensg_to_gene, enst_to_isoname, ensp_to_gene and isoname_to_gene files have been prepared")

def IsoLenTab(fa_file, results):
    """
    Prepare a table that provides gene and length information for each isoform
    """

    # Initialize Lists
    isos = []
    genes = []
    lens = []

    # Parse fafsa file to append isoform, gene and length information 
    for line in open(fa_file):
        if line.startswith('>'):
            isos.append(line.split('|')[4].split('""')[0])
            genes.append(line.split('|')[5].split('""')[0])
            lens.append(line.split('|')[6].split('""')[0])

    # Export Data as a DataFrame and a tsv file
    data = {'isoform': isos, 'gene': genes, 'length': [int(x_len) for x_len in lens]}
    df = pd.DataFrame(data)
    df.to_csv(results.isoname_lens, sep='\t', index=False)
    print("The isoform length table has been prepared.")
    return df

def GeneLenTab(IsolenFile, results):
    """ 
    Prepare a table that provides the min, max and average length of a gene 
    """

    cut = IsolenFile[['gene', 'length']]

    # Map genes to lengths and calc average, min and max. Round mean to nearest tenth. Reset indices. 
    length= cut.groupby(['gene']).length.agg(['mean', 'min', 'max'])
    length['mean'] = length['mean'].round(decimals = 1)
    length = length.reset_index(level=['gene'])

    # Change column names and save the table 
    length.columns =['gene', 'avg_len', 'min_len', 'max_len']
    length.to_csv(results.gene_lens, sep="\t", index=False)
    print('Prepared the gene length statistics table')

    
def protein_coding_genes(results):
    df_gtf = gtfparse.read_gtf(results.gtf_file)
    df_gtf = df_gtf[df_gtf['feature'] == 'gene']
    df_gtf = df_gtf[df_gtf['gene_type'] == 'protein_coding']
    protein_coding_genes = df_gtf['gene_name'].unique()
    with open(results.pc_genes, 'w') as filehandle:
        for gene in protein_coding_genes:
            filehandle.write(f"{gene}\n")

def main():
    
    # If results folder does not exist, make it 
    # TODO - maybe remove
    rdir = './dump'
    if not os.path.exists(rdir):
        os.mkdir(rdir)


    # Command line arguments
    parser = argparse.ArgumentParser(description='Proccess ORF related file locations')
    parser.add_argument('--gtf','-g',action='store', dest= 'gtf_file',help='Gencode GTF input file location')
    parser.add_argument('--fa','-fa',action='store', dest= 'fa_file',help='Gencode Fafsa input file location')
    parser.add_argument('--ensg_gene', '-oeg', action='store', dest='ensg_gene', help='ensg to gene output file location') 
    parser.add_argument('--enst_isoname', '-oei', action='store', dest='enst_isoname', help='enst to isoname output file location')
    parser.add_argument('--gene_ensp', '-oge', action='store', dest='gene_ensp', help='Gene to ensp output file location')
    parser.add_argument('--gene_isoname', '-ogi', action='store', dest='gene_isoname', help="Gene to isoname output file location") 
    parser.add_argument('--isoname_lens', '-oigl', action='store', dest='isoname_lens', help='Isoname length table output location')
    parser.add_argument('--gene_lens', '-ogls', action='store', dest='gene_lens', help='Gene Length statistics output location')
    parser.add_argument('--protein_coding_genes', '-pcg', action='store', dest='pc_genes', help='Protein Coding genes output location')
    results = parser.parse_args()

    # Make ensg -> gene, enst -> isoname, ensp -> gene and isoname -> gene mapping files 
    GenMap(results.gtf_file, results)

    # Prepare Gene Isoform Length table 
    df = IsoLenTab(results.fa_file, results)

    # Prepare Gene Length Table 
    GeneLenTab(df, results)
    
    protein_coding_genes(results)
    
    

if __name__ == "__main__":
    main()


