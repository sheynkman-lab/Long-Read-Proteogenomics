#!/usr/bin/env python3
import pandas as pd 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--protein_classification',action='store',dest='protein_classification')
parser.add_argument('--best_orf',action='store',dest='best_orf')
parser.add_argument('--refined_meta',action='store',dest='refined_meta')
parser.add_argument('--name',action='store',dest='name')
parser.add_argument('--dest_dir',action='store',dest='dest_dir')
args = parser.parse_args()

pclass = pd.read_table(args.protein_classification)
best_orf = pd.read_table(args.best_orf)
refined = pd.read_table(args.refined_meta)
pclass = (
    pclass
        .merge(best_orf[['pb_acc','orf_calling_confidence','has_stop_codon']], left_on='pb',right_on='pb_acc', how='inner')
        .merge(refined[['base_acc','gene','CPM',]], left_on='pb_acc',right_on='base_acc', how='inner')
        .drop(columns=['pb_acc','base_acc'])
)
pclass.to_csv(f'{args.dest_dir}{args.name}.protein_classification_w_meta.tsv', sep='\t',index=False)