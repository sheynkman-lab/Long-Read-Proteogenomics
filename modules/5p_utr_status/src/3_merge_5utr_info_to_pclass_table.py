#!/usr/bin/env python3
# merge 5utr info into the main protein classification table

#%%

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--name',action='store',dest='name')
parser.add_argument('--utr_info',action='store',dest='utr_info')
parser.add_argument('--sqanti_protein_classification',action='store',dest='sqanti_protein_classification')
parser.add_argument('--odir',action='store',dest='odir')
args = parser.parse_args()

pclass = pd.read_table(args.sqanti_protein_classification)

utr_info = pd.read_table(args.utr_info)
utr_info = utr_info[['pb', 'utr_exon_status', 'utr_cat']]

pclass = pd.merge(pclass, utr_info, how='left', on='pb')

pclass.to_csv(os.path.join(args.odir,f'{args.name}.sqanti_protein_classification_w_5utr_info.tsv'), sep='\t', index=None)
