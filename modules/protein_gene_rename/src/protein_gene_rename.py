#!/usr/bin/env python3

import argparse
import pandas as pd 
from Bio import SeqIO

def rename_fasta_genes(sample_protein_fasta, pb_gene, name):
    renamed_fasta = []
    for record in SeqIO.parse(sample_protein_fasta, 'fasta'):
        acc = record.description.split('|')[1]
        if acc in pb_gene.keys():
            record.description = f'{record.id} GN={pb_gene[acc]}'
            renamed_fasta.append(record)
    SeqIO.write(renamed_fasta, f'{name}.protein_refined.fasta','fasta')

def rename_gtf_genes(sample_gtf, pb_gene, name):
    with open(sample_gtf, 'r') as ifile, open(f'{name}_with_cds_refined.gtf','w') as cds_ofile:
        for line in ifile:
            seqname,source,feature,start,end,score,strand,phase,attributes=line.split('\t')
            gene_info,pb_acc,cpm = attributes.split("|")
            if pb_acc in pb_gene.keys():
                gene = pb_gene[pb_acc]
                gene_info = f'gene_id "{gene}"; transcript_id "{gene}'
                new_attributes='|'.join([gene_info,pb_acc,cpm])
                new_line='\t'.join([seqname,source,feature,start,end,score,strand,phase, new_attributes])
                cds_ofile.write(new_line)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_gtf',action='store',dest='sample_gtf')
    parser.add_argument('--sample_protein_fasta',action='store',dest='sample_fasta')
    parser.add_argument('--pb_protein_genes',action='store',dest='pb_protein_genes')
    parser.add_argument('--name',action='store',dest='name')
    args = parser.parse_args()

    pb_gene = pd.read_table(args.pb_protein_genes)
    pb_gene = pd.Series(pb_gene.pr_gene.values,index=pb_gene.pb).to_dict()

    rename_fasta_genes(args.sample_fasta, pb_gene, args.name)
    rename_gtf_genes(args.sample_gtf, pb_gene, args.name)




if __name__=="__main__":
    main()