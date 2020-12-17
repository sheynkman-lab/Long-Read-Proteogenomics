#!/usr/bin/env python3

from gtfparse import read_gtf
from collections import defaultdict
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
import logging

def orf_mapping(orf_coord, gencode, sample_gtf, orf_seq):
    def get_num_upstream_atgs(row):
        orf_start = int(row['orf_start'])
        acc = row['pb_acc']
        if acc in orf_seq.keys():
            seq = orf_seq[acc] 
        else:
            return np.inf
        upstream_seq = seq[0:orf_start-1] # sequence up to the predicted ATG
        num_atgs = upstream_seq.count('ATG')
        return num_atgs
    
    exons = sample_gtf[sample_gtf['feature'] == 'exon'].copy()
    exons['exon_length'] = abs(exons['end'] - exons['start']) + 1
    exons.rename(columns = {'start' : 'exon_start', 'end': 'exon_end'}, inplace = True)
    start_codons = gencode[gencode['feature'] == 'start_codon']
    start_codons = start_codons[['seqname','transcript_id','strand',  'start', 'end']].copy()
    start_codons.rename(columns = {'start' : 'start_codon_start', 'end': 'start_codon_end'}, inplace = True)
    
    logging.info(f"exons\n{exons}")
    logging.info(f"start_codons\n{start_codons}")
    
    plus = plus_mapping(exons, orf_coord, start_codons[start_codons['strand'] == '+'])
    logging.info(f"Plus mapping complete\n{plus.head()}")
    
    minus = minus_mapping(exons, orf_coord, start_codons[start_codons['strand'] == '-'])
    all_cds = pd.concat([plus, minus])
    all_cds['upstream_atgs'] = all_cds.apply(get_num_upstream_atgs, axis=1)
    return all_cds

def plus_mapping(exons, orf_coord, start_codons):
    """
    Map plus strand
    """
    def compare_start_plus(row, start_codons):
        start = int(row['cds_start'])
        match = start_codons[start_codons['start_codon_start'] == start]
        match = match[match['seqname'].str.strip() == row['seqname'].strip()]
        if len(match) > 0:
            return list(match['transcript_id'])
        return None

    plus_exons = exons[exons['strand'] == '+'].copy()
    if len(plus_exons) == 0:
        logging.info("no plus exons")
        return plus_exons
    plus_exons['current_size'] = plus_exons.sort_values(by = ['transcript_id', 'exon_start']).groupby('transcript_id')['exon_length'].cumsum()
    plus_exons['prior_size'] = plus_exons['current_size'] - plus_exons['exon_length']

    plus_comb = pd.merge(orf_coord, plus_exons, left_on = 'pb_acc', right_on = 'transcript_id', how = 'inner')
    if len(plus_comb) ==0:
        logging.info("No matches found between cpat orfs and sample exons")
        return plus_comb
    logging.info(f"Combined\n{plus_comb.head()}")
    plus_comb = plus_comb[(plus_comb['prior_size'] <= plus_comb['orf_start']) &  (plus_comb['orf_start'] <= plus_comb['current_size'])]
    if len(plus_comb) ==0:
        logging.info("No matching start CDS exon")

    plus_comb['start_diff'] = plus_comb['orf_start'] - plus_comb['prior_size']
    plus_comb['cds_start'] = plus_comb['exon_start'] + plus_comb['start_diff'] - 1
    plus_comb['gencode_atg'] = plus_comb.apply(lambda row : compare_start_plus(row, start_codons), axis = 1)
    plus_comb.drop(columns=['exon_length', 'current_size', 'prior_size', 'start_diff'], inplace = True)
    return plus_comb
    

def minus_mapping(exons, orf_coord, start_codons):
    """
    Map minus strand
    """
    def compare_start_minus(row, start_codons):
        start = int(row['cds_start'])
        match = start_codons[(start_codons['start_codon_end'] == start) ]
        match = match[match['seqname'].str.strip() == row['seqname'].strip()]
        if len(match) > 0:
            return list(match['transcript_id'])
        return None
    
    minus_exons = exons[exons['strand'] == '-'].copy()
    
    minus_exons['current_size'] = minus_exons.sort_values(by = ['transcript_id', 'exon_start'], ascending=[True,False]).groupby('transcript_id')['exon_length'].cumsum()
    minus_exons['prior_size'] = minus_exons['current_size'] - minus_exons['exon_length']
    minus_comb = pd.merge(orf_coord, minus_exons, left_on = 'pb_acc', right_on = 'transcript_id', how = 'inner')
    if len(minus_comb) == 0:
        return minus_comb
    minus_comb = minus_comb[minus_comb['orf_start'].between(minus_comb['prior_size'],minus_comb['current_size'])]
    minus_comb['start_diff'] = minus_comb['orf_start'] - minus_comb['prior_size']
    minus_comb['cds_start'] = minus_comb['exon_end'] - minus_comb['start_diff'] + 1
    minus_comb['gencode_atg'] = minus_comb.apply(lambda row : compare_start_minus(row, start_codons), axis = 1)
    minus_comb.drop(columns=['exon_length', 'current_size', 'prior_size', 'start_diff'], inplace = True)
    return minus_comb

def read_orf(filename):
    """
    Reads the ORF file and formats the column names
    Keep only predictions with score higher than protein-coding threshold

    Parameters
    ---------
    filename : str
        location of orf coordinate file

    Returns
    --------
    orf: pandas DataFrame
    """
    orf = pd.read_csv(filename, sep = '\t')
    orf[['pb_acc', 'misc', 'orf']] = orf['ID'].str.split('_', expand=True)
    orf['pb_acc'] = orf['pb_acc'].str.split('|').str[0]
    orf = orf.drop(labels=['ID', 'misc'], axis=1)
    orf.columns = ['len', 'orf_strand', 'orf_frame', 'orf_start', 'orf_end', 'orf_len', 'fickett', 'hexamer', 'coding_score', 'pb_acc', 'orf_rank',]
    logging.info(f"ORF file read \n{orf.head()}")
    return orf

def orf_calling(orf, num_orfs_per_accession = 1):
    """
    Choose 'best' ORF by examining number of upstream ATG's, match to Gencode Transcript start, and ORF codings
    """
    def call_orf(group):
        score_threshold = 0.364
        def calling_confidence(row):
            highscore = 0.9
            if row['atg_rank'] == 1 and row['score_rank'] == 1:
                return 'Clear Best ORF'
            elif row['coding_score'] > score_threshold and row['atg_rank'] == 1:
                return 'Nonsense Mediated Decay'
            elif row['coding_score'] <= score_threshold:
                return 'Low Quality ORF'
            return 'Plausable ORF'
        
        group['atg_rank'] = group['upstream_atgs'].rank(ascending=True)
        group['score_rank'] = group['coding_score'].rank(ascending=False)
        group['orf_calling_confidence'] = group.apply(lambda row : calling_confidence(row), axis = 1)
        
        with_gencode = group.dropna(subset=['gencode_atg'])
        if len(with_gencode) >=1:
            group = with_gencode
        
#         good_score = group[group['coding_score'] >= score_threshold]
#         if(len(good_score) >= 1):
#             group = good_score
        atg_shift = 10       # how much to shift sigmoid for atg score
        group['atg_score'] = group['upstream_atgs'].apply(lambda x : 1 - 1/( 1+ np.exp(-x+atg_shift))
        group['atg_score'] = group['upstream_atgs'].apply(lambda x : 1/x  if x > 1 else 0.99)
#         group['gencode_score'] = group['gencode_atg'].apply(lambda x : 0 if x == '' else 0.8)
        group['orf_score'] = group.apply(lambda row: 1 - (1-row['coding_score']*0.95)*(1-row['atg_score']), axis = 1)
#         group = group.sort_values(by='atg_rank').reset_index(drop=True)
#         if group.loc[0,'atg_rank'] == group.loc[0,'score_rank']:
#             return group.head(1)
        
        group = group.sort_values(by='orf_score', ascending=False).reset_index(drop=True)
        return group.head(num_orfs_per_accession)
        
    called_orf = orf.groupby('pb_acc').apply(call_orf).reset_index(drop=True)
    return called_orf



def orf_calling_v1(orf):
    """
    Choose 'best' ORF by examining number of upstream ATG's, match to Gencode Transcript start, and ORF codings
    """
    def call_orf(group):
        def calling_confidence(row):
            if row['atg_rank'] == 1 and row['score_rank'] == 1:
                return 'Clear Best ORF'
            elif row['coding_score'] > highscore and row['atg_rank'] == 1:
                return 'Nonsense Mediated Decay'
            elif row['coding_score'] < score_threshold:
                return 'Low Quality ORF'
            return 'None'

        score_threshold = 0.364

        good_score = group[group['coding_score'] >= score_threshold]
        if(len(good_score) >= 1):
            group = good_score
            
        highscore = 0.9

        group['orf_calling_confidence'] = 'None'
        group['atg_rank'] = group['upstream_atgs'].rank(ascending=True)
        group['score_rank'] = group['coding_score'].rank(ascending=False)

        group['atg_score'] = group['upstream_atgs'].apply(lambda x : 1/x  if x > 1 else 0.99)
        group['gencode_score'] = group['gencode_atg'].apply(lambda x : 0 if x == '' else 0.8)
        group['orf_score'] = group.apply(lambda row: 1 - (1-row['coding_score']*0.99)*(1-row['atg_score'])*(1-row['gencode_score']), axis = 1)

        

        group['orf_calling_confidence'] = group.apply(lambda row : calling_confidence(row), axis = 1)
        group = group.sort_values(by='atg_rank').reset_index(drop=True)
        if group.loc[0,'atg_rank'] == group.loc[0,'score_rank']:
            return group.head(1)
        
        group = group.sort_values(by='orf_score', ascending=False).reset_index(drop=True)
        return group.head(1)
        
    called_orf = orf.groupby('pb_acc').apply(call_orf).reset_index(drop=True)
    return called_orf
    
    
def main():
    parser = argparse.ArgumentParser(description='Proccess ORF related file locations')
    parser.add_argument('--orf_coord', '-oc',action='store', dest= 'orf_coord',help='ORF coordinate input file location')
    parser.add_argument('--gencode_gtf','-g',action='store', dest= 'gencode_gtf',help='gencode coordinate input file location')
    parser.add_argument('--sample_gtf','-sg',action='store', dest= 'sample_gtf',help='Sample GTF input file location')
    parser.add_argument('--pb_gene','-pg',action='store', dest= 'pb_gene',help='PB Accession/Gencode id mapping input file location')
    parser.add_argument('--classification','-c',action='store', dest= 'classification',help='sample classification input file location')
    parser.add_argument('--sample_fasta','-sf',action='store', dest= 'sample_fasta',help='Sample FASTA input file location')
    parser.add_argument('--output','-o',action='store', dest= 'output',help='Output file location')
    results = parser.parse_args()

    orf_coord = read_orf(results.orf_coord)
    gencode = read_gtf(results.gencode_gtf)
    sample_gtf = read_gtf(results.sample_gtf)
    pb_gene = pd.read_csv(results.pb_gene, sep = '\t')
    classification = pd.read_csv(results.classification, sep = '\t')

    orf_seq= defaultdict() # pb_acc -> orf_seq
    for rec in SeqIO.parse(results.sample_fasta, 'fasta'):
        pb_id = rec.id.split('|')[0] 
        orf_seq[pb_id] = str(rec.seq)


    all_orfs = orf_mapping(orf_coord, gencode, sample_gtf, orf_seq)
    orfs = orf_calling(all_orfs, 100)

    classification = classification[['isoform', 'FL']]
    total = classification['FL'].sum()
    classification['CPM'] = classification['FL'] / total * 1000000

    orfs = pd.merge(orfs, pb_gene, left_on = 'pb_acc', right_on='isoform', how = 'left')
    orfs = pd.merge(orfs, classification, on = 'isoform', how = 'left')
    orfs = orfs.drop(columns = ['isoform'])
    orfs.to_csv(results.output, index = False, sep = "\t")


if __name__ == "__main__":
    main()
