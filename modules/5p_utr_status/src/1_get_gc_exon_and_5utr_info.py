#!/usr/bin/env python3
# distinguish trunc versus novel for tx nic/nnc
# read in gencode gtf, only exons from protein-coding genes, and make a merged bed
# then find out if protruding 5' end is within intron or exons to determine if
# an nic/nnc is subset (and creating an ntrunc) or protruding (and creating a novel atg start)

#%%

import subprocess
from collections import Counter
from collections import defaultdict 
import os
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gencode_gtf',action='store',dest='gencode_gtf')
parser.add_argument('--odir',action='store',dest='odir')
args = parser.parse_args()



odir = args.odir
if not os.path.exists(odir):
    os.mkdir(odir)

# prepare files and data structures for later comparison of pacbio 5utr to gencode

# read in isonames containing a cds (coding)
pc_isonames = set() # isonames with a cds
biotypes = []
for line in open(args.gencode_gtf):
    if line.startswith('#'): continue
    wds = line.split('\t')
    if wds[2] == 'CDS':
        isoname = wds[8].split('transcript_name "')[1].split('"')[0]
        pc_isonames.add(isoname)

# make a bed file with cds isonames
gencode_cds_bed_fpath = os.path.join(odir,'gencode_exons_for_cds_containing_ensts.bed')
with open(gencode_cds_bed_fpath, 'w') as ofile: 
    for line in open(args.gencode_gtf):
        if line.startswith('#'): continue
        wds = line.split('\t')
        if wds[2] == 'exon':
            isoname = wds[8].split('transcript_name "')[1].split('"')[0]
            if isoname in pc_isonames:
                chrom, start, end = wds[0], wds[3], wds[4]
                ofile.write('\t'.join([chrom, start, end]) + '\n')

# make gencode merged bed
# this will be used for determining if 5' end of pacbio transcripts are protruding into intronic regions
cat_cmd = f'''cat {gencode_cds_bed_fpath} '''
grep_cmd = '''| bedtools sort | bedtools merge | uniq | grep '^chr' | awk '{print $1 "\t" ($2 - 1) "\t"  $3}' '''
out_cmd = f'''> {gencode_cds_bed_fpath}_merged.bed'''
cmd = cat_cmd + grep_cmd + out_cmd
process = subprocess.run(cmd, shell=True, capture_output=True) 

# load in gc exon coords
gc_coords = defaultdict(lambda: defaultdict(lambda: [None, None, []])) # genename -> [<strand>, <exon chain string>, <exon coords as ints>]
# example - GAPDH -> ['+', '50-100_150-200_250-300', [[50, 100], [150, 200], [250, 300]]]
for line in open(args.gencode_gtf):
    if line.startswith('#'): continue
    wds = line.split('\t')
    if wds[2] == 'exon':
        isoname = wds[8].split('transcript_name "')[1].split('"')[0]
        genename = wds[8].split('gene_name "')[1].split('"')[0]
        if isoname in pc_isonames:
            chrom, start, end, strand = wds[0], int(wds[3]), int(wds[4]), wds[6]
            gc_coords[genename][isoname][0] = strand
            gc_coords[genename][isoname][2].append([start, end])

# make exon chains
gc_chains = defaultdict(lambda: defaultdict(lambda: [])) # gene -> isoname -> exon chain string
for gene, iso_dict in gc_coords.items():
    for isoname, [strand, _, exon_coords] in iso_dict.items():
        exon_coords = sorted(exon_coords)
        coord_str = '_'.join([str(s) + '-' + str(e) for s, e in exon_coords])
        gc_coords[gene][isoname][1] = coord_str
        gc_chains[gene][isoname] = coord_str


# write out exon chain strings
with open(os.path.join(odir, 'gc_exon_chain_strings_for_cds_containing_transcripts.tsv'), 'w') as ofile:
    for gene in gc_chains:
        for isoname in gc_chains[gene]:
            chain_string = gc_chains[gene][isoname]
            ofile.write(f'{gene}\t{isoname}\t{chain_string}\n')



