#!/usr/bin/env python3

# make GTF file that includes the ORF regions (as CDS features)
# input - pacbio gtf ('jurkat_corrected.gtf'), orf calls ('jurkat_refine_orf_calls.tsv')
# output - pacbio gtf with added "cds" features (orfs)

# %%

import pandas as pd
from collections import defaultdict
import copy
import argparse
import gtfparse



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

def ceiling_cpm(cpm):
    # gtf top score is 1000
    if cpm > 1000:
        return 1000
    else:
        return cpm




def make_pacbio_cds_gtf(sample_gtf, agg_orfs, refined_orfs, pb_gene, output_cds):
    """
    sample_gtf : filename
    agg_orfs : filename
    refined_orfs : filename
    pb_gene : filename
    output_cds : filename
    """
    # import gtf, only exon info.
    # only move forward with representative pb isoform (for same-protein groups)
    gtf = gtfparse.read_gtf(sample_gtf)

    gtf = gtf[['seqname', 'feature', 'start', 'end', 'strand', 'gene_id']]
    gtf = gtf[gtf['feature'] == 'exon']
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
    # only move forward with "base accession" (representative pb)
    repr_pbs = pd.read_table(agg_orfs)['base_acc'].to_list()
    gtf = gtf[gtf.acc.isin(repr_pbs)]


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
    ranges = pd.read_table(refined_orfs)[['pb_acc', 'orf_start', 'orf_end', 'CPM']]

    # read in pb to genename
    pb_gene = pd.read_table(pb_gene)
    pb_gene = pd.Series(pb_gene.gene.values, index=pb_gene.pb_acc).to_dict()


    with open(output_cds, 'w') as ofile:
        for i, row in ranges.iterrows():
            acc, orf_start, orf_end, cpm = row
            if acc in pbs:
                gene = pb_gene[acc]
                # only continue with pb isoforms that align to a genetic locus
                if gene == '-': continue
                infos = pbs[acc]
                chr, strand, coords, blens, cblens = infos
                # NOTE - uncomment to do chr22-oly test 
                # if chr != 'chr22': continue
                i1, delta1 = get_first_block_index(orf_start, cblens)
                i2, delta2 = get_first_block_index(orf_end, cblens)
                if strand == '+':
                    orf_coords = make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)
                elif strand == '-':
                    orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords)
                # write out the coordinates
                acc_w_gene_w_cpm = gene + '|' + acc + '|' + str(cpm)
                out_acc = 'gene_id "{}"; transcript_id "{}";'.format(gene, acc_w_gene_w_cpm)
                for [start, end] in coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'exon', str(start), str(end), '.', strand, '.', out_acc]) + '\n')
                for [start, end] in orf_coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'CDS', str(start), str(end), '.', strand, '.', out_acc]) + '\n')



def main():
    parser = argparse.ArgumentParser("IO file locations for make pacbio cds gtf")
    parser.add_argument("--sample_gtf", action="store", dest = "sample_gtf")
    parser.add_argument("--agg_orfs", action="store", dest = "agg_orfs")
    parser.add_argument("--refined_orfs", action="store", dest="refined_orfs")
    parser.add_argument("--pb_gene", action="store", dest="pb_gene")
    parser.add_argument("--output_cds", action="store", dest="output_cds")
    results = parser.parse_args()
    make_pacbio_cds_gtf(results.sample_gtf, results.agg_orfs, results.refined_orfs, results.pb_gene, results.output_cds)
    
if __name__ == "__main__":
    main()

