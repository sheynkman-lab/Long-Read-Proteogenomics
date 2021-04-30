#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import pandas as pd 

parser = argparse.ArgumentParser()
parser.add_argument('--name',action='store',dest='name')
parser.add_argument('--collapsed_fasta',action='store',dest='collapsed_fasta')
parser.add_argument('--classification',action='store',dest='classification')
args = parser.parse_args()
collapsed_accs = set()
for record in SeqIO.parse(args.collapsed_fasta,'fasta'):
    collapsed_accs.add(record.id)

classification = pd.read_table(args.classification)
classification = classification[classification['isoform'].isin(collapsed_accs)]
classification.to_csv(f'{args.name}_classification.5degfilter.tsv',sep='\t',index=False)