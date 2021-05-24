#!/usr/bin/env python3

# derive a same-protein-sequence gencode file


import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import argparse

def make_gencode_database(gencode_fasta, output_fasta, output_cluster, protein_coding_genes):
    # make gencode clusters of same-protein-sequence entries
    # evaluate protein sequence similarity within the same gene
    gc = defaultdict(lambda: defaultdict(list)) # prot_seq -> gene -> [isonames]
    for rec in SeqIO.parse(gencode_fasta, 'fasta'):
        id_split = rec.id.split('|')
        gene = id_split[6]
        isoname = id_split[5]
        if gene in protein_coding_genes:
            gc[rec.seq][gene].append(isoname)

    # write a fasta, with same-protein sequences clustered
    # for each group of same-protein-sequence isonames, choose the first alphanum
    with open(output_fasta, 'w') as ofile, open(output_cluster, 'w') as ofile2:
        ofile2.write('repr_isoname\tsame_prot_seq_isonames\n')
        for seq in gc:
            for gene, isonames in gc[seq].items():
                isonames = sorted(isonames)
                repr_isoname = isonames[0]
                ofile.write('>gc|{}| GN={}\n{}\n'.format(repr_isoname, gene, seq))
                ofile2.write(repr_isoname + '\t' + ','.join(isonames[1:]) + '\n')
    


def main():
    parser = argparse.ArgumentParser("Makes Gencode Clusters of same protein sequence entries")
    parser.add_argument('--gencode_fasta', '-if', action='store', dest='gencode_fasta', help = 'gencode fasta file to group')
    parser.add_argument('--protein_coding_genes', '-pc', action='store',dest='protein_coding_genes', help='protein coding genes to keep in fasta file')
    parser.add_argument('--output_fasta', '-of', action='store',dest='output_fasta', help = 'output file location of fasta file')
    parser.add_argument('--output_cluster', '-oc', action='store', dest='output_cluster', help = 'output cluster tsv file location')
    results = parser.parse_args()

    gencode_fasta = results.gencode_fasta
    output_fasta = results.output_fasta
    output_cluster = results.output_cluster
    with open(results.protein_coding_genes, 'r') as file:
        protein_coding_genes = file.read().splitlines()
    make_gencode_database(gencode_fasta, output_fasta, output_cluster, protein_coding_genes)

if __name__ == "__main__":
    main()

#%%
# %%
