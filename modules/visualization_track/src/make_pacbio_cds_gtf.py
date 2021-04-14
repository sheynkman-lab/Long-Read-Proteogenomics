#!/usr/bin/env python3

# make GTF file that includes the ORF regions (as CDS features)
# input - pacbio gtf ('jurkat.collapsed.gtf'), orf calls ('jurkat_refine_orf_calls.tsv')
# output - pacbio gtf with added "cds" features (orfs)

# %%

import pandas as pd
import numpy as np
from collections import defaultdict
import copy
import argparse
import gtfparse
import logging

logger = logging.getLogger('cds_logger')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('make_pacbio_cds.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
ch.setFormatter(formatter)
fh.setFormatter(formatter)
# add the handlers to logger
logger.addHandler(ch)
logger.addHandler(fh)

def string_to_boolean(string):
    """
    Converts string to boolean

    Parameters
    ----------
    string :str
    input string

    Returns
    ----------
    result : bool
    output boolean
    """
    if isinstance(string, bool):
        return str
    if string.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif string.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_first_block_index(orf_coord, cblens, pblens):
    # get the index corresponding to the first block containing the orf start
    # return index, and the dela (spacing upstream of end)
    for i, cblen in enumerate(cblens):
        if orf_coord <= cblen:
            delta = cblen - orf_coord
            return i, delta
    logger.warning(f"ORF COORDINATE IS NOT FOUND WITHIN BLOCKS")
    return i, 0

def make_cds_coords_positive_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][0] = orf_coords[0][1] - delta1
    orf_coords[-1][1] = orf_coords[-1][1] - delta2
    return orf_coords

def make_cds_coords_negative_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][1] = orf_coords[0][0] + delta1 
    orf_coords[-1][0] = orf_coords[-1][0] + delta2 
    return orf_coords

# def get_first_block_index(orf_coord, cblens, pblens):
#     """
#     Finds the index of the block where orf_coord location is and the relative 
#     difference (delta) within the block the coordinate is
    
#     Args:
#         orf_coord (int): orf coordinate location
#         cblens (numpy 1D array) : cumulative block lengths
#         pblens (numpy 1D array ) : prior cumulative block lengths
    
#     Returns:
#         i (int): index of orf location
#         delta (int): difference within block orf_coord starts
    
#     """
#     # get the index corresponding to the first block containing the orf start
#     # return index, and the dela (spacing upstream of end)
#     # logger.info(f"orf_coord {orf_coord} \ncblens {cblens}")
#     for i, (cblen, pblen) in enumerate(zip(cblens, pblens)):
#         if orf_coord <= cblen:
#             delta = orf_coord - pblen
#             return i, delta
#     logging.WARNING(f"ORF COORDINATE IS NOT FOUND WITHIN BLOCKS")
#     return i, cblen - pblen

# def make_cds_coords_positive_strand(i1, delta1, i2, delta2, coords):
#     """
#     Makes the CDS coordinates for the positive strand

#     Args:
#         i1 (int): index of first CDS exon
#         delta1 (int) : offset from start of first exon to make CDS
#         i2 (int) : index of last CDS exon
#         delta2 (int) : delta of end from start - length of last CDS exon
#         coords [[int,int]] : coordinates of exons

#     Returns:
#         [[int,int]]: CDS coordinates
    
#     """
#     logger.info(f"\n+\n{i1}\t{delta1}\n{i2}\t{delta2}\n{coords}")
#     orf_coords = copy.deepcopy(coords)
#     orf_coords = orf_coords[i1: i2+1]
#     # trim ends to orf start/end
#     orf_coords[0][0] = orf_coords[0][0] + delta1
#     orf_coords[-1][1] = orf_coords[-1][0] + delta2
#     logger.info(f"\n{orf_coords}")
#     return orf_coords

# def make_cds_coords_negative_strand(i1, delta1, i2, delta2, coords):
    # """Makes the CDS coordinates for the negative strand

    # Args:
    #     i1 (int): index of start CDS exon
    #     delta1 (int): offset from end of first exon to make CDS
    #     i2 ([type]): index of last CDS exon
    #     delta2 ([type]): offset from last CDS exon end as start location
    #     coords ([[int,int]]):  coordinates of exons

    # Returns:
    #     [[int,int]]: CDS coordinates
    # """

    # logger.info(f"\n-\n{i1}\t{delta1}\n{i2}\t{delta2}\n{coords}")
    # orf_coords = copy.deepcopy(coords)
    # orf_coords = orf_coords[i1: i2+1]
    # # trim ends to orf start/end
    # orf_coords[0][1] = orf_coords[0][1] - delta1 + 2
    # orf_coords[-1][0] = orf_coords[-1][1] - delta2 + 2
    # logger.info(f"\n{orf_coords}")
    # return orf_coords

def ceiling_cpm(cpm, ceiling = 1000):
    """Sets CPM to have ceiling

    Args:
        cpm (float): CPM of isoform
        ceiling (int, optional): Maximum. Defaults to 1000.

    Returns:
        float: new cpm constrained to ceiling
    """
    # gtf top score is 1000
    if cpm > ceiling:
        return ceiling
    else:
        return cpm

#%%
def get_min_and_max_coords_from_exon_chain(coords):
    """Gets the minumum and maximum coordinates from exon chain

    Args:
        coords ([[int, int]]): [[start,end]] exon coordinate chain

    Returns:    
        (int, int): start and end coordinates of entire chain
    """
   
    
    min_coord = min(map(min, coords))
    max_coord = max(map(max, coords))
    return min_coord, max_coord
    

#%%
def make_pacbio_cds_gtf(sample_gtf, refined_orfs, called_orfs, pb_gene, name, include_transcript):
    """Makes PacBio CDS and saves file with CDS

    Args:
        sample_gtf (filename): sample_gtf file 
        refined_orfs (filename): aggregate_orf info. from Refined DB
        called_orfs (filename): orf calls from ORF_Calling
        pb_gene (filename): PacBio gene name cross reference
        name (string): name of sample
        include_transcript (bool): whether to include transcript in saved gtf file 
    """
    refined_db = pd.read_table(refined_orfs)
    representative_accessions = refined_db['base_acc'].to_list()
    # import gtf, only exon info.
    # only move forward with representative pb isoform (for same-protein groups)
    gtf = gtfparse.read_gtf(sample_gtf)

    gtf = gtf[['seqname', 'feature', 'start', 'end', 'strand', 'transcript_id']]
    gtf = gtf[gtf['feature'] == 'exon']
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
    # only move forward with "base accession" (representative pb)

    gtf = gtf[gtf.acc.isin(representative_accessions)]

    # pb coords into dict
    pbs = defaultdict(lambda: ['chr', 'strand', [], [], [],[]]) # pb -> [chr, strand, [start, end], [block lengths],[cum. block lengths], [prior cumulative block lengths]]
    # PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
    for i, row in gtf.iterrows():
        chr, feat, start, end, strand, acc = row
        pbs[acc][0] = chr
        pbs[acc][1] = strand
        pbs[acc][2].append([int(start), int(end)])
    # sort all coords, calc blocks
    for acc, infos in pbs.items():
        strand = infos[1]
        if strand == '+':
            infos[2] = sorted(infos[2])
        elif strand == '-':
            infos[2] = sorted(infos[2], reverse=True)
        infos[3] = np.array([end-start+1 for [start, end] in infos[2]])
        infos[4] = np.cumsum(infos[3])
        infos[5] = infos[4] - infos[3]


    # read in the ranges of orf on pb transcripts
    ranges = pd.read_table(called_orfs)[['pb_acc', 'orf_start', 'orf_end']]
    ranges = pd.merge(
        ranges, refined_db[['base_acc', 'CPM']], 
        left_on='pb_acc',
        right_on='base_acc',
        how='inner')
    ranges = ranges[['pb_acc', 'orf_start', 'orf_end', 'CPM']]
    # ranges = ranges[ranges['pb_acc'].isin(representative_accessions)]

    # read in pb to genename
    pb_gene = pd.read_table(pb_gene)
    pb_gene = pd.Series(pb_gene.gene.values, index=pb_gene.pb_acc).to_dict()


    with open(f"{name}_with_cds.gtf", 'w') as ofile:
        for i, row in ranges.iterrows():
            acc, orf_start, orf_end, cpm = row
            cpm = round(cpm)
            # remove stop exon
            orf_end = orf_end - 3
            if acc in pbs:
                if acc in pb_gene.keys():
                    gene = pb_gene[acc]
                else: 
                    gene = '-'
                # only continue with pb isoforms that align to a genetic locus
                if gene == '-': continue
                
                infos = pbs[acc]
                chr, strand, coords, blens, cblens, pblens = infos
                # NOTE - uncomment to do chr22-oly test 
                # if chr != 'chr22': continue

                logger.info(f"\n{acc}\t{strand}\n{orf_start} - {orf_end}\n{blens}\n{cblens}")
                
                
                i1, delta1 = get_first_block_index(orf_start, cblens, pblens)
                i2, delta2 = get_first_block_index(orf_end, cblens, pblens)
                if strand == '+':
                    orf_coords = make_cds_coords_positive_strand(i1, delta1, i2, delta2, coords)
                elif strand == '-':
                    orf_coords = make_cds_coords_negative_strand(i1, delta1, i2, delta2, coords)
                # write out the coordinates
                acc_w_gene_w_cpm = gene + '|' + acc + '|' + str(cpm)
                out_acc = f'gene_id "{gene}"; transcript_id "{acc_w_gene_w_cpm}";'

                if include_transcript:
                    tstart, tend = get_min_and_max_coords_from_exon_chain(coords)
                    ofile.write('\t'.join([chr, 'hg38_canon', 'transcript', str(tstart), str(tend), '.', strand, '.', out_acc]) + '\n')
                for [start, end] in coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'exon', str(start), str(end), '.', strand, '.', out_acc]) + '\n')
                for [start, end] in orf_coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'CDS', str(start), str(end), '.', strand, '.', out_acc]) + '\n')



def main():
    
    parser = argparse.ArgumentParser("IO file locations for make pacbio cds gtf")
    parser.add_argument("--name", action="store", dest="name", help="name of sample - used for output file name")
    parser.add_argument("--sample_gtf", action="store", dest = "sample_gtf", help="sample gtf, from sqanti3")
    parser.add_argument("--called_orfs", action="store", dest = "called_orfs", help="agg orf tsv file, from refined DB")
    parser.add_argument("--refined_database", action="store", dest="refined_database",help="best orfs per PB accession, from ORFCalling" )
    parser.add_argument("--pb_gene", action="store", dest="pb_gene", help="filename of pb-gene reference file")
    parser.add_argument("--include_transcript", action="store", dest="include_transcript", help="yes/no - include transcript in output file")
    results = parser.parse_args()
    include_transcript = string_to_boolean(results.include_transcript)
    make_pacbio_cds_gtf(results.sample_gtf, results.refined_database, results.called_orfs, results.pb_gene, results.name, include_transcript)


if __name__ == "__main__":
    main()
    
"""
python src/make_pacbio_cds_gtf.py \
--name test \
--sample_gtf   test_data/test.collapsed.gff \
--agg_orfs     test_data/test.orf.aggregated.tsv \
--refined_orfs test_data/test.orf_calls.tsv \
--pb_gene      test_data/test.pb_gene.tsv \
--include_transcript no
"""


