#!/usr/bin/env python3
# make a six frame translation database for from the pacbio transcripts

from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import argparse

# read in pacbio transcripts

def return_all_orfs(translated_seq):
    # return all contiguous polypeptide sequences 7 AA or higher
    seqs = translated_seq.split('*')
    filtered_seqs = []
    for seq in seqs:
        if len(seq) >= 7:
            filtered_seqs.append(str(seq))
    return filtered_seqs

def make_pacbio6fm_gene_grouped(iso_annot, ensg_gene, sample_fasta, output_fasta):
    """
    Runs PacBio 6 Frame Gene Grouping and saves results

    Parameters
    -----------
    iso_annot : filename
        sample classification file from SQANTI3
    ensg_gene : filename
        ensg_gene file from ReferenceTables module
    sample_fasta : filename :
        sample corrected fasta file
    output_fasta : filename
        output filename of grouped pacbio results
    """
    # get associated gene for each pb acc, for final write-out (below)
    # pb_gene is pb_acc -> gene dictionary
    df = pd.read_table(iso_annot)[['isoform', 'associated_gene']]
    df2 = pd.read_table(ensg_gene, header=None)
    df2.columns = ['associated_gene', 'gene']
    df3 = pd.merge(df, df2, on='associated_gene', how='left').fillna('NOVEL')[['isoform', 'gene']]
    pb_gene = pd.Series(df3.gene.values, index=df3.isoform).to_dict()

    gene_seqs = defaultdict(lambda: set()) # gene -> pacbio sequences as list

    for rec in SeqIO.parse(sample_fasta, 'fasta'):
        pb_id = rec.id.split('|')[0] 
        gene = pb_gene[pb_id]
        F1 = rec.seq.translate()
        F2 = rec.seq[1:].translate()
        F3 = rec.seq[2:].translate()
        R1 = rec.seq.reverse_complement().translate()
        R2 = rec.seq.reverse_complement()[1:].translate()
        R3 = rec.seq.reverse_complement()[2:].translate()
        translations = [F1, F2, F3, R1, R2, R3]
        for tr in translations:
            orfs = set(return_all_orfs(tr))
            gene_seqs[gene].update(orfs)


    # write out fasta file in which entries represent each gene and the
    # pacbio-derived protein "space"
    with open(output_fasta, 'w') as ofile:
        for gene, orfs in gene_seqs.items():
            ofile.write('>' + gene + '\n' + '-'.join(orfs) + '\n')


def main():
    parser = argparse.ArgumentParser(description='PacBIO 6 frame gene grouping')
    parser.add_argument('--iso_annot', action='store', dest='iso_annot', help = 'SQANTI3 jurkat classification file location')
    parser.add_argument('--ensg_gene', action='store', dest='ensg_gene',help='ENSG_GENE file location from ReferenceTables')
    parser.add_argument('--sample_fasta', action='store', dest='sample_fasta', help='Sample corrected fasta file')
    parser.add_argument('--output_fasta', '--output', action='store', dest='output_fasta', help = 'output: grouped fasta file')
    results = parser.parse_args()
    make_pacbio6fm_gene_grouped(results.iso_annot, results.ensg_gene, results.sample_fasta, results.output_fasta)

if __name__ == "__main__":
    main()