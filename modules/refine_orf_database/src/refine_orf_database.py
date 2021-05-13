#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio import Seq
from collections import defaultdict
import argparse
import logging
import os


def get_accession_seqs(seqs):
    logging.info('getting accesssion sequences...')
    pb_seqs = defaultdict() # pb_acc -> transcript_seq
    redundant_accs = []
    for entry in seqs:
        seq = str(entry.seq)
        pb_acc = entry.id.split('|')[0]
        # skip duplicates
        if pb_acc in pb_seqs:
            redundant_accs.append(pb_acc)
        else:
            pb_seqs[pb_acc] = seq
    return pb_seqs,redundant_accs
    
    
def combine_by_sequence(orfs, pb_seqs):
    logging.info('combining by sequence...')
    orfs = orfs[['pb_acc', 'orf_start', 'orf_end', 'orf_len']]
    # extract, translate, and aggregate protein sequences
    pb_pseqs = defaultdict(lambda: list()) # protein_seq -> list of pb acc
    for index, row in orfs.iterrows():
        pb_acc, start, end, olen = list(row)
        seq = pb_seqs[pb_acc]
        orf_seq = seq[start-1:end]
        prot_seq = Seq.translate(orf_seq, to_stop=True)
        pb_pseqs[prot_seq].append(pb_acc)
    return pb_pseqs


def order_pb_acc_numerically(accs):
    logging.info('Ordering PB Accession Numerically...')
    # order pb accessions by numbers
    accs_numerical = []
    for acc in accs:
        pb, gene_idx, iso_idx = acc.split('.')
        gene_idx = int(gene_idx)
        iso_idx = int(iso_idx)
        accs_numerical.append([gene_idx, iso_idx])
    accs_numerical.sort()
    num_sorted_accs = ['PB.' + str(g) + '.' + str(i) for g, i in accs_numerical]
    return num_sorted_accs



    

def aggregate_results(pacbio, orfs):
    logging.info('Aggregating Results')
    def get_total(accessions, orf_dict):
        total = 0
        for acc in accessions:
            total += orf_dict[acc]
        return total
    
    pacbio['accessions'] = pacbio['pb_accs'].str.split('|')
    pacbio['base_acc'] = pacbio['accessions'].apply(lambda x : x[0])
    
    fl_dict = pd.Series(orfs.FL.values,index=orfs.pb_acc).to_dict()
    cpm_dict = pd.Series(orfs.CPM.values,index=orfs.pb_acc).to_dict()
    
    pacbio['FL'] = pacbio['accessions'].apply(lambda accs : get_total(accs, fl_dict))
    pacbio['CPM'] = pacbio['accessions'].apply(lambda accs : get_total(accs, cpm_dict))
    return pacbio

def filter_orf_scores(orfs, cutoff):
    logging.info(f'Filtering ORF coding_score to be above {cutoff}')
    """
    Filter ORFS so only orfs above a cutoff value are used.
    Per CPAT a score >= 0.364 is considered a protein-coding score
    orfs : pandas DataFrame
        called ORFs
    cutoff : float
        remove all orfs that have coding_score < cutoff
    """
    orfs = orfs[orfs['coding_score'] >= cutoff]
    return orfs

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

    
def main():
    parser = argparse.ArgumentParser(description='Proccess ORF related file locations')
    parser.add_argument('--name', action='store', dest='name', help='sample name')
    parser.add_argument('-io', '--orfs',action='store', dest= 'orfs',help='ORF coordinate input file location')
    parser.add_argument('-if', '--pb_fasta', action = 'store', dest='pb_fasta', help='PacBio fasta sequence input file location')
    parser.add_argument('-cut', '--coding_score_cutoff', dest = 'cutoff', type = float, default = 0.0, help='CPAT coding score cutoff. remove all below')
    results = parser.parse_args()
    name = results.name
    logging.info('Reading Fasta File...')
    seqs = SeqIO.parse(open(results.pb_fasta), 'fasta')
    logging.info('Reading ORFS...')
    orfs = pd.read_table(results.orfs)
    
    # Filter ORFS based on score and whether protein coding
    orfs = filter_orf_scores(orfs, results.cutoff)

    # only keep orfs that have a stop codon
    orfs = orfs.query('has_stop_codon')
    
    pb_seqs,redundant_accs = get_accession_seqs(seqs)
    pb_pseqs = combine_by_sequence(orfs, pb_seqs)
    
    # Save combined results

    with open(f'{name}_combined.tsv', 'w') as ofile, open(f'{name}_combined.fasta', 'w') as ofile2:
        ofile.write('protein_sequence\tpb_accs\n')
        for seq, accs in pb_pseqs.items():
            accs_sorted = order_pb_acc_numerically(accs)
            accs_str = '|'.join(accs_sorted)
            ofile.write(seq + '\t' + accs_str + '\n')
            ofile2.write('>' + accs_str + '\n' + seq + '\n')
    
    
    pacbio = pd.read_csv(f'{name}_combined.tsv', sep = '\t')
    seqs = SeqIO.parse(open(f'{name}_combined.fasta'), 'fasta')
    os.remove(f'{name}_combined.tsv')
    os.remove(f'{name}_combined.fasta')

    
    pacbio = aggregate_results(pacbio, orfs)
    orfs = orfs[['pb_acc', 'coding_score', 'orf_score', 'orf_calling_confidence', 'upstream_atgs', 'gene']]
    pacbio = pd.merge(pacbio, orfs, how = 'inner', left_on = 'base_acc', right_on = 'pb_acc')
    logging.info('Writing refined database fasta results...')
    pb_gene = pd.Series(orfs.gene.values,index=orfs.pb_acc).to_dict()
    base_map = pd.Series(pacbio.base_acc.values,index=pacbio.pb_accs).to_dict()
    with open(f'{name}_orf_refined.fasta', 'w') as ofile:
        for entry in seqs:
            seq = str(entry.seq)
            pb_acc = entry.id
            base_acc = base_map[pb_acc]
            gene = pb_gene[base_acc]
            ofile.write(f'>pb|{base_acc}|fullname GN={gene}\n{seq}\n')
    
    logging.info('Writing refined datavase tsv results...')
    pacbio = pacbio[['pb_accs', 'base_acc', 'coding_score','orf_calling_confidence','upstream_atgs','orf_score','gene','FL', 'CPM']]
    pacbio.to_csv(f'{name}_orf_refined.tsv', sep = '\t', index = False) 
    logging.info('Refine Database Complete')
    logging.info('************************')
    
    
if __name__ == '__main__':
    main()