#!/usr/bin/env python3

# derive a same-protein-sequence gencode file

# %%

import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import argparse




def make_gencode_database(gencode_fasta, ouptut_fasta, output_cluster):
    # make gencode clusters of same-protein-sequence entries
    # evaluate protein sequence similarity within the same gene
    gc = defaultdict(lambda: defaultdict(list)) # prot_seq -> gene -> [isonames]
    for rec in SeqIO.parse(gencode_fasta, 'fasta'):
        gene = rec.id.split('|')[6]
        isoname = rec.id.split('|')[5]
        gc[rec.seq][gene].append(isoname)

    # write a fasta, with same-protein sequences clustered
    # for each group of same-protein-sequence isonames, choose the one 
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
    parser.add_argument('--output_fasta', '-of', action='store',dest='output_fasta', help = 'output file location of fasta file')
    parser.add_argument('--output_cluster', '-oc', action='store', dest='output_cluster', 'output cluster tsv file location')
    gencode_fasta = results.gencode_fasta
    output_fasta = results.output_fasta
    output_cluster = results.output_cluster
    make_gencode_database(gencode_fasta, output_fasta, output_cluster)

