#!/usr/bin/env python3

#%%
import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bed',action='store',dest='bed')
parser.add_argument('--name',action='store',dest='name')
args = parser.parse_args()
name = args.name

bed_names = ['chrom','chromStart','chromStop','acc_full','score','strand','thickStart','thickEnd','itemRGB','blockCount','blockSizes','blockStarts']
bed = pd.read_table(args.bed, names=bed_names)
bed['acc_full'] = bed['acc_full'].str.split('(').str[0]

with open(f'{name}_peptides.bed12', 'w') as ofile:
        ofile.write(f'track name={name}_peptides\n')
        bed.to_csv(ofile, sep='\t', index=None, header=None)


bed['rgb'] = '0,51,0'
filter_names = ['chrom','chromStart','chromStop','acc_full','score','strand','thickStart','thickEnd','rgb','blockCount','blockSizes','blockStarts']
bed = bed[filter_names]
with open(f'{name}_shaded_peptides.bed12', 'w') as ofile:
    ofile.write(f'track name="{name} PacBio Peptides" itemRgb=On\n')
    bed.to_csv(ofile, sep='\t', index=None, header=None)
