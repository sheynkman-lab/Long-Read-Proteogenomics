# make GTF file including the ORF regions (as CDS features)
# before writing out the GTF file, cluster same-orf-structure entries

import pandas as pd
from collections import defaultdict
import copy

# import gtf, only exon info.
gtf = pd.read_table('../../a_SQANTI3_out/jurkat_corrected.gtf', header=None)
gtf['acc'] = gtf[8].apply(lambda x: x.split('gene_id "')[1].split('"')[0])
gtf = gtf[[0, 2, 3, 4, 6, 'acc']]
gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
gtf = gtf.loc[gtf['feat']=='exon']


# pb coords into dict
pbs = defaultdict(lambda: ['chr', 'strand', [], [], []]) # pb -> [chr, strand, [start, end], [block lengths], [cum. block lengths]]
# PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
for i, row in gtf.iterrows():
    chr, feat, start, end, strand, acc = row
    pbs[acc][0] = chr
    pbs[acc][1] = strand
    pbs[acc][2].append([int(start), int(end)])



# sort all coords, calc blocks
def make_cumulative_blens(blocks):
    cblocks = []
    cbl = 0 # cumulative block length
    for b in blocks:
        cbl += b
        cblocks.append(cbl)
    return cblocks

for acc, infos in pbs.items():
    strand = infos[1]
    if strand == '+':
        infos[2] = sorted(infos[2])
    elif strand == '-':
        infos[2] = sorted(infos[2], reverse=True)
    infos[3] = [end-start+1 for [start, end] in infos[2]]
    infos[4] = make_cumulative_blens(infos[3])



# read in the ranges of orf on pb transcripts
ranges = pd.read_table('a_jurkat_cpat_orf_ranges.tsv')



# trim coords to make orf-range entries
def get_first_block_index(orf_coord, cblens):
    # get the index corresponding to the first block containing the orf start
    # return index, and the dela (spacing upstream of end)
    for i, cblen in enumerate(cblens):
        if orf_coord <= cblen:
            delta = cblen - orf_coord
            return i, delta

def make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][0] = orf_coords[0][1] - delta1
    orf_coords[-1][1] = orf_coords[-1][1] - delta2
    return orf_coords

def make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][1] = orf_coords[0][0] + delta1
    orf_coords[-1][0] = orf_coords[-1][0] + delta2
    return orf_coords

with open('a_gtf_w_cds.gtf', 'w') as ofile:
    for i, row in ranges.iterrows():
        acc, orf_start, orf_end = row
        if acc in pbs:
            infos = pbs[acc]
            chr, strand, coords, blens, cblens = infos
            i1, delta1 = get_first_block_index(orf_start, cblens)
            # print(acc, orf_start, orf_end)
            i2, delta2 = get_first_block_index(orf_end, cblens)
            if strand == '+':
                orf_coords = make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)
            elif strand == '-':
                orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords)
            # write out the coordinates
            for [start, end] in coords:
                ofile.write('\t'.join([chr, 'hg38_canon', 'exon', str(start), str(end), '.', strand, '.', acc]) + '\n')
            for [start, end] in orf_coords:
                ofile.write('\t'.join([chr, 'hg38_canon', 'CDS', str(start), str(end), '.', strand, '.', acc]) + '\n')
