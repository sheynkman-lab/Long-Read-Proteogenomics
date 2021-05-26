#!/usr/bin/env python3
# based on identified peptides, make a bed file

# %%

import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import copy
import argparse
import logging


def most_frequent(List):
	counter = 0
	element = List[0]
	for i in List:
		curr_frequency = List.count(i)
		if(curr_frequency > counter):
			counter = curr_frequency
			element = i
	return element


# sort all coords, calc blocks
def make_cumulative_blens(blocks):
    cblocks = []
    cbl = 0 # cumulative block length
    for b in blocks:
        cbl += b
        cblocks.append(cbl)
    return cblocks


def read_sample_gtf(gtf_file):
    gtf = pd.read_table(gtf_file, skiprows=1, header=None)
    gtf = gtf[[0, 2, 3, 4, 6, 8]]
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
    gtf = gtf.loc[gtf['feat']=='CDS']
    gtf['acc'] = gtf['acc'].str.split('|').str[1]
    return gtf

def read_reference_gtf(gtf_file):
    gtf = pd.read_table(gtf_file, skiprows=5, header=None)
    gtf = gtf[[0, 2, 3, 4, 6, 8]]
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
    gtf = gtf.loc[gtf['feat']=='CDS']

    gtf['acc'] = gtf['acc'].str.split('transcript_name "').str[1]
    gtf['acc'] = gtf['acc'].str.split('";').str[0]
    return gtf

def process_gtf(gtf):
    # CDS coords into dict
    pbs = defaultdict(lambda: ['chr', 'strand', [], [], []]) # pb -> [chr, strand, [start, end], [block lengths], [cum. block lengths]]
    # PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
    for i, row in gtf.iterrows():
        chr, feat, start, end, strand, acc = row
        pbs[acc][0] = chr
        pbs[acc][1] = strand
        pbs[acc][2].append([int(start), int(end)])
    
    for acc, infos in pbs.items():
        strand = infos[1]
        if strand == '+':
            infos[2] = sorted(infos[2])
        elif strand == '-':
            infos[2] = sorted(infos[2], reverse=True)
        infos[3] = [end-start+1 for [start, end] in infos[2]]
        infos[4] = make_cumulative_blens(infos[3])
    return pbs


def process_psmtsv(psmtsv_file, gene_map):

    # read in peptide data, filter for target + 1% fdr
    pep = pd.read_table(psmtsv_file)
    # pep = pep.iloc[:, [12, 24, 37, 48]]
    pep = pep[['Base Sequence','Protein Accession','Decoy/Contaminant/Target','QValue', 'Previous Amino Acid', 'Next Amino Acid']]
    pep.columns = ['pep', 'pb_acc', 'targ', 'qval', 'prev_aa','next_aa']
    # TODO add back in proper filter
    # pep = pep[pep.targ == 'T']
    pep = pep[(pep.targ =='T') & (pep.qval <= 0.01)]
    # for now, all accessions
    # TODO - decide how to deal with multiple PBs mapped to peptide
    pep['pb_acc'] = pep['pb_acc'].str.split('|')
    pep = pep[['pep', 'pb_acc', 'prev_aa','next_aa']]
    pep = pep.explode('pb_acc')

    # add in gene
    pep['gene'] = pep['pb_acc'].map(gene_map)
    return pep


# #1 determine number of proteins to which each peptide maps (specificity score)
# #2 for each peptide, determine the first and last index when mapped to the pb_acc CDS

# load pb_acc to protein sequence mapping
# have ~112K protein sequences (from 175K pb transcript sequences)

def read_fasta(fasta_file):
    seqs = defaultdict() # pb_acc -> protein_sequence
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        pb_acc = rec.id.split('|')[1]
        seqs[pb_acc] = str(rec.seq)
    return seqs



# using the table of pb_acc, peptide, start, end, determine the coordinates


def get_first_block_index(orf_coord, cblens):
    # get the index corresponding to the first block containing the orf start
    # return index, and the dela (spacing upstream of end)
    for i, cblen in enumerate(cblens):
        if orf_coord <= cblen:
            delta = cblen - orf_coord
            return i, delta
    return i, 0

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


def write_peptide_gtf(name, pep_ranges, pbs, gene_pb, seqs):
    # Note - code below was modified from script that found CDS orf range wtihin
    # a full-length RNA isoforms, so orf-range is now peptide range 
    with open(f'{name}_peptides.gtf', 'w') as ofile:
        # write UCSC track header
        # remove for conversion to bed12 (genePred complains)
        # ofile.write('track name=peptide color=0,0,0\n')
        for i, row in pep_ranges.iterrows():
            
            pep_seq, pb_acc, prev_aa, next_aa,gene, num_prots, pep_start, pep_end = row
            # convert from protein (AA) to CDS (nt) coords
            pep_start = pep_start * 3 - 2
            pep_end = pep_end * 3
            if pb_acc in pbs:
                infos = pbs[pb_acc]
                chr, strand, coords, blens, cblens = infos
                i1, delta1 = get_first_block_index(pep_start, cblens)
                i2, delta2 = get_first_block_index(pep_end, cblens)
                if strand == '+':
                    orf_coords = make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)
                elif strand == '-':
                    orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords)
                # write out the coordinates
                prev_aa = most_frequent(prev_aa.split('|'))
                next_aa = most_frequent(next_aa.split('|'))
                if chr in ['chrX','chrY']:
                    gene = f"{gene}_{chr}"
                acc_id= f"{prev_aa}.{pep_seq}.{next_aa}({gene})"
                pep_acc = f'gene_id "{acc_id}"; transcript_id "{acc_id}";'
                for [start, end] in orf_coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'CDS', str(start), str(end), '.', strand,
                                '.', pep_acc]) + '\n')

def make_peptide_gtf(name , pbs, pb_gene, pep, seqs):
    def find_start_end_pep_index(row):
        pep, pb_acc = row[['pep', 'pb_acc']]
        # print(f"{pb_acc}\t{pep}")
        seq = seqs[pb_acc]
        start = seq.find(pep)
        if start == -1:
            return
        start_idx = start + 1
        end_idx = start_idx + len(pep) - 1
        row['start'] = start_idx
        row['end'] = end_idx
        return row
    
    def get_number_of_mapping_proteins(row):
        pep = row['pep']
        gene = row['gene']
        if pd.isnull(gene):
            return 0
        pb_accs = gene_pb[gene] 
        num_prots = 0
        for pb_acc in pb_accs:
            if pb_acc not in seqs:
                # non-base pacbios, ignore
                continue
            if pep in seqs[pb_acc]:
                num_prots += 1
        return num_prots
    # gene to pacbio ditionary
    gene_pb = {}
    for pb_acc, gene in pb_gene.items():
        gene_pb[gene] = gene_pb.get(gene, []) + [pb_acc]
    
    pep['num_prots'] = pep.apply(get_number_of_mapping_proteins, axis = 1)
    pep = pep.apply(find_start_end_pep_index, axis=1)
    pep = pep.drop_duplicates()
    write_peptide_gtf(name, pep, pbs, gene_pb, seqs)


def main():
    parser = argparse.ArgumentParser("IO file locations for make peptide gtf file")
    parser.add_argument("--name", action="store", dest="name", help="name of sample - used for output file name")
    parser.add_argument("--sample_gtf", action="store", dest = "sample_gtf", help="sample gtf with cds. from make_pacbio_cds_gtf")
    parser.add_argument("--reference_gtf",action="store",dest="reference_gtf")
    parser.add_argument("--peptides", action="store", dest="peptides", help = "peptides file location. from MetaMorpheus")
    parser.add_argument("--pb_gene", action="store", dest="pb_gene", help ="PB-Gene reference")
    parser.add_argument("--gene_isoname", action="store",dest="gene_isoname")
    parser.add_argument("--refined_fasta", action="store", dest="refined_fasta", help = "refined fasta file. from refined db generation")
    results = parser.parse_args()
    sample_gtf = read_sample_gtf(results.sample_gtf)
    sample_pbs = process_gtf(sample_gtf)

    reference_gtf = read_reference_gtf(results.reference_gtf)
    reference_pbs = process_gtf(reference_gtf)
    pbs = {**sample_pbs, **reference_pbs}
    pb_gene = pd.read_table(results.pb_gene)
    pb_gene = pd.Series(pb_gene.pr_gene.values, index=pb_gene.pb).to_dict()
    gene_isoname = pd.read_table(results.gene_isoname, names=['gene','isoname'])
    gene_isoname = pd.Series(gene_isoname.gene.values, index=gene_isoname.isoname).to_dict()
    gene_map = {**pb_gene, **gene_isoname}
    pep = process_psmtsv(results.peptides, gene_map)
    seqs = read_fasta(results.refined_fasta)
    make_peptide_gtf(results.name, pbs, gene_map, pep, seqs)

#%%


if __name__ == "__main__":
    main()
# %%
# """
# python src/make_peptide_gtf_file.py \
# --name test \
# --sample_gtf test_data/test_with_cds.gtf \
# --peptides test_data/test.peptides.psmtsv \
# --pb_gene test_data/test.pb_gene.tsv \
# --refined_fasta test_data/test.orf.aggregated.fasta
# """

# #%%
# pb_gene = pd.read_table("/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/huvec/pb_gene.tsv")
# pb_gene = pd.Series(pb_gene.gene.values, index=pb_gene.pb_acc).to_dict()
# pep = process_psmtsv("/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/huvec/huvec_psms.tsv", pb_gene)
# # %%
# gtf, pbs = process_gtf('/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/make_cds/results/huvec/pacbio_cds/huvec_with_cds.gtf')
# # %%
# seqs = read_fasta("/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Visualization/data/huvec/huvec_orf_aggregated.fasta")
# # %%
#%%
# make_peptide_gtf('test', gtf, pbs, pb_gene, pep, seqs)
# %%
