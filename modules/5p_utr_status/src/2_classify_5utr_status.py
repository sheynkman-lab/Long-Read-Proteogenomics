#!/usr/bin/env python3
# classify 5utr status for pb proteins 
# classification focuses on monoexonic and multiexonic 5utrs

#%%

from collections import defaultdict
from collections import Counter
import re
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--gencode_exons_bed',action='store',dest='gencode_exons_bed')
parser.add_argument('--gencode_exons_chain',action='store',dest='gencode_exons_chain')
parser.add_argument('--sample_cds_gtf', action='store',dest='sample_cds_gtf')
parser.add_argument('--odir',action='store',dest='odir')
args = parser.parse_args()




## read in gencode merged exon ranges (ranges for exons coming from ensts that have cds)
gc_exons = defaultdict(lambda: []) # chr -> list of exon ranges
# read in gencode bed file
for line in open(args.gencode_exons_bed):
    chrom, start, end = line.split('\t')
    start, end = int(start), int(end)
    gc_exons[chrom].append([start, end])

## read in gencode coords
gc_chains = defaultdict(lambda: []) # gene -> list of exon chains (str)
for line in open(args.gencode_exons_chain):
    gene, isoname, exon_chain_str = line.strip().split('\t')
    gc_chains[gene].append(exon_chain_str)


## read in pb info into defaultdict

pb_info = defaultdict(lambda: [None, None, None, -1, -1, -1, None])
# pb -> [gene, chrom, strand, tss, nterm, num_5utr_exons, 5utr junction string]

#region - reading in pb info defauldict info

# read in pb info - temporary
pbs = defaultdict(lambda: [None, None, None, [], []])
# pb -> [gene, chrom, strand, exons, cdss]
for line in open(args.sample_cds_gtf):
    gene = line.split('gene_id "')[1].split('"')[0]
    pb = line.split('|')[1]
    wds = line.split('\t')
    chrom, strand, start, end = wds[0], wds[6], int(wds[3]), int(wds[4])
    pbs[pb][0] = gene 
    pbs[pb][1] = chrom
    pbs[pb][2] = strand
    if wds[2] == 'exon':
        pbs[pb][3].append([start, end])
    elif wds[2] == 'CDS':
        pbs[pb][4].append([start, end])


## start helper functions

def get_the_tss_coord(exons, strand):
    # assumes sorted exon coords
    if strand == '+':
        return exons[0][0]
    else:
        return exons[-1][1]

def get_the_nterm_coord(cdss, strand):
    # assumes sorted cds coords
    if strand == '+':
        return cdss[0][0]
    else:
        return cdss[-1][1]

def derive_5utr_coords(nterm, exons, strand):
    coords_5utr = []
    if strand == '+':
        for s, e in exons:
            if s <= nterm <= e:
                coords_5utr.append([s, nterm])
                break
            else:
                coords_5utr.append([s, e])
    else: # strand == '-'
        for s, e in sorted(exons, reverse=True):
            if s <= nterm <= e:
                coords_5utr.append([nterm, e])
                break
            else:
                coords_5utr.append([s, e])
    return sorted(coords_5utr)

def derive_junction_chain_from_exon_chain(exon_chain_str):
    # assume format is start-end_start-end
    if '_' not in exon_chain_str:
        # only one exon, return nothing
        return ''
    else:
        internal_coords = exon_chain_str.split('-')[1:-1]
        internal_coords_str = '-'.join(internal_coords)
        return internal_coords_str

# read in pb info (see pb_info defaultdict above)
for pb, [gene, chrom, strand, exons, cdss] in pbs.items():
    exons = sorted(exons)
    cdss = sorted(cdss)
    tss = get_the_tss_coord(exons, strand)
    nterm = get_the_tss_coord(cdss, strand)
    coords_5utr = derive_5utr_coords(nterm, exons, strand)
    num_5utr_exons = len(coords_5utr)
    coords_5utr_str = '_'.join([str(s) + '-' + str(e) for s, e in coords_5utr])
    junc_str_5utr = derive_junction_chain_from_exon_chain(coords_5utr_str)
    pb_info[pb] = [gene, chrom, strand, tss, nterm, num_5utr_exons, junc_str_5utr]


# determine 5utr properties, write out to table


# helper functions

def get_5end_coord_at_pb_junction_match(strand, junc_chain_5utr, gc_chain):
    # get 5' end of the 5' most exon that mapped to the pb 5utr junction chain
    if strand == '+':
        prefix = gc_chain.split(junc_chain_5utr)[0]
        linked_5end_coord = int(re.split('_|-', prefix)[-2])
        return linked_5end_coord
    else: # neg strand
        prefix = gc_chain.split(junc_chain_5utr)[1]
        linked_5end_coord = int(re.split('_|-', prefix)[1])
        return linked_5end_coord

def determine_if_5end_of_junction_chain_is_protruding_or_subset(strand, tss, junc_chain_5utr, gc_chain):
    # status would be protrudes or is_subset
    # require protruding by 10 or more nt to be a valid protrusion
    gc_5end = get_5end_coord_at_pb_junction_match(strand, junc_chain_5utr, gc_chain)
    if strand == '+':
        if (tss + 9) < gc_5end:
            return 'protruding'
        else:
            return 'subset'
    else: # neg strand
        if (tss - 9) > gc_5end:
            return 'protruding'
        else:
            return 'subset'

def get_5utr_junc_chain_status(gene, strand, tss, junc_chain_5utr, gc_chains):
    # determine 5'utr category, if 5utr has one or more junctions
    # first determine if junction chain is known or novel
    junc_stat = 'novel'
    tss_stat = ''
    for gc_chain in gc_chains[gene]:
        # gc_chain is the junction chain for the whole gc transcript (e.g., '100-150_200-250_300-375')
        if junc_chain_5utr in gc_chain:
            tss_stat = determine_if_5end_of_junction_chain_is_protruding_or_subset(strand, tss, junc_chain_5utr, gc_chain)
            if tss_stat == 'subset':
                junc_stat = 'perfect_subset'
            else: # protruding
                if junc_stat != 'perfect_subset':
                    junc_stat = 'known_protruding'
    return junc_stat

def get_exon_status(num_5utr_exons):
    exon_status = 'monoexonic'
    if num_5utr_exons > 1:
        exon_status = 'multiexonic'
    return exon_status

def get_5utr_mono_exon_status(gene, strand, tss, gc_ranges):
    # note that if protrusion is not equal to or greater than 10 nt
    # that the tss is considered to be "within" the gc exons
    for (start, end) in gc_ranges:
        if start <= tss <= end:
            # tss resides within gc exons
            return 'within'
    # at this point, tss is in intronic region
    # determine the protruding distance
    protruding_distance = 999999999 # chr1 len is 250M
    if strand == '+':
        for (start, end) in gc_ranges:
            if start > tss:
                cur_dist = start - tss
                if cur_dist < protruding_distance:
                    protruding_distance = cur_dist
    else: # strand == '-'
        for (start, end) in gc_ranges:
            if tss > end:
                cur_dist = tss - end
                if cur_dist < protruding_distance:
                    protruding_distance = cur_dist
    if protruding_distance >= 10:
        return 'protruding'
    else:
        return 'within'

def is_pb_nterm_within_gc_exons(pb_start, gc_ranges):
    for (start, end) in gc_ranges:
        if start <= pb_start <= end:
            return True
    return False

def determine_final_utr_cat(junc_cat):
    if 'protruding' in junc_cat or 'novel' in junc_cat:
        return 'unique'
    else: # monoexon within or multiexonic perfect_subset
        return 'subset'


with open(os.path.join(args.odir, 'pb_5utr_categories.tsv'), 'w') as ofile:
    ofile.write('pb\tgene\tnum_5utr_exons\tutr_exon_status\ttss_in_gc_exons\tjunc_cat\tutr_cat\n')
    lens = []
    for pb, [gene, chrom, strand, tss, nterm, num_5utr_exons, junc_chain_5utr] in pb_info.items():
        gc_ranges = gc_exons[chrom]
        junc_cat = '-'
        exon_status = get_exon_status(num_5utr_exons)
        if exon_status == 'monoexonic':
            junc_cat = get_5utr_mono_exon_status(gene, strand, tss, gc_ranges)
        else: # exon_status == 'multiexonic':
            junc_cat = get_5utr_junc_chain_status(gene, strand, tss, junc_chain_5utr, gc_chains)
        is_tss_within_gc_exons = is_pb_nterm_within_gc_exons(tss, gc_ranges) * 1
        utr_cat = determine_final_utr_cat(junc_cat)
        odata = [pb, gene, num_5utr_exons, exon_status, is_tss_within_gc_exons, junc_cat, utr_cat]
        ofile.write('\t'.join([str(x) for x in odata]) + '\n')
        
#%%

import matplotlib.pyplot as plt
lens2 = sorted(lens)[0:1000]
plt.hist(lens, bins=range(50))







