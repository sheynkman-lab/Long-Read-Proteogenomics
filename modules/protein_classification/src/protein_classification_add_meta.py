#!/usr/bin/env python3
import pandas as pd 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--protein_classification',action='store',dest='protein_classification')
parser.add_argument('--best_orf',action='store',dest='best_orf')
parser.add_argument('--refined_meta',action='store',dest='refined_meta')
parser.add_argument('--ensg_gene',action='store',dest='ensg_gene')
parser.add_argument('--name',action='store',dest='name')
parser.add_argument('--dest_dir',action='store',dest='dest_dir')
args = parser.parse_args()

pclass = pd.read_table(args.protein_classification)
best_orf = pd.read_table(args.best_orf)
refined = pd.read_table(args.refined_meta)
ensg_gene = pd.read_table(args.ensg_gene, names=['ENSG','Gene'])
ensg_gene = pd.Series(ensg_gene.Gene.values,index=ensg_gene.ENSG).to_dict()

pclass = (
    pclass
        .merge(best_orf[['pb_acc','orf_calling_confidence','has_stop_codon']], left_on='pb',right_on='pb_acc', how='inner')
        .merge(refined[['base_acc','gene','CPM',]], left_on='pb_acc',right_on='base_acc', how='inner')
        
)

def get_gene_names(ensgs, gene_dict):
    gene_names = set()
    for e in ensgs:
        if e in gene_dict:
            gene_names.add(gene_dict[e])
        else:
            gene_names.add(e)
        
    return list(gene_names)

pclass['pr_genes'] = pclass['pr_gene'].str.split(',')
pclass['pr_gene_names'] = pclass['pr_genes'].apply(lambda ensgs: get_gene_names(ensgs,ensg_gene))

def get_matching_gene(row):
    tx_gene = row['gene']
    for pr_gene in row['pr_gene_names']:
        if pr_gene==tx_gene:
            return pr_gene
    return row['pr_gene_names'][0]
pclass['pr_gene'] = pclass.apply(get_matching_gene, axis = 1) 
pclass['tx_gene'] = pclass['gene']
pclass = pclass.drop(columns=['pb_acc','base_acc', 'pr_genes', 'pr_gene_names', 'gene'])

pclass[['pb','tx_gene','pr_gene']].to_csv(f'{args.dest_dir}{args.name}_genes.tsv',sep='\t',index=False)
pclass.to_csv(f'{args.dest_dir}{args.name}.protein_classification_w_meta.tsv', sep='\t',index=False)

