#!/usr/bin/env python3

# temp scripts for different parts of the pipline

import pandas as pd
from Bio import SeqIO
from Bio import Seq
from collections import defaultdict
import argparse




def get_accession_seqs(seqs):
    pb_seqs = defaultdict() # pb_acc -> transcript_seq
    redundant_accs = []
    for entry in seqs:
        seq = str(entry.seq)
        pb_acc = entry.id
        # skip duplicates
        if pb_acc in pb_seqs:
            redundant_accs.append(pb_acc)
        else:
            pb_seqs[pb_acc] = seq
    return pb_seqs,redundant_accs
    
    
def combine_by_sequence(orfs, pb_seqs):
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


def process_args():
    parser = argparse.ArgumentParser(description='Proccess ORF related file locations')
    parser.add_argument('-io', '--orfs',action='store', dest= 'orfs',help='ORF coordinate input file location')
    parser.add_argument('-if', '--pb_fasta', action = 'store', dest='pb_fasta', help='PacBio fasta sequence input file location')
    parser.add_argument('-ipc', '--protien_coding', action = 'store', dest='protein_coding_file', help='Gencode protien coding genes input file location')
    parser.add_argument('-or', '--redundant', action='store', dest='red', help = 'Output redundant accession file location')
    parser.add_argument('-oct', '--combined_tsv', action='store', dest='combined_tsv', help = 'Output combined tsv file location')
    parser.add_argument('-ocf', '--combined_fasta', action='store', dest='combined_fasta', help = 'Output combined fasta file location')
    parser.add_argument('-ot', '--agg_tsv', action='store', dest='agg_tsv', help = 'Output aggregated tsv file location')
    parser.add_argument('-of', '--agg_fasta', action='store', dest='agg_fasta', help = 'Output aggregated fasta file location')
    parser.add_argument('-pc', '--protein_coding_only', dest= 'protein_coding_only', default = 'no', help ='Keep only protein coding genes', )
    parser.add_argument('-cut', '--coding_score_cutoff', dest = 'coding_score_cutoff', 
                        type = 'float', default = 0.0, help='CPAT coding score cutoff. remove all below')
    results = parser.parse_args()
    return results
    

def aggregate_results(pacbio, orfs):
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
    """
    Filter ORFS so only orfs above a cutoff value are used.
    Per CPAT a score >= 0.364 is considered a protein-coding score
    orfs : pandass DataFrame
        called ORFs
    cutoff : float
        remove all orfs that have coding_score < cutoff
    """
    orfs = orfs[orfs['coding_score'] >= cutoff]
    return orfs

def filter_orf_protein_coding(orfs, protein_coding_filename):
    """
    Filter ORFs to only contain genes that are known to be protein-coding
    per Gencode
    
    Parameters
    ----------
    orfs : pandass DataFrame
        called ORFs
    protein_coding_filename : filename
        file of protien-coding genes. text file seperated by lines
    """
    with open(protein_coding_filename, 'r') as file:
        protein_coding_genes = file.read().splitlines()
    orfs = orfs[orfs['gene'].isin(protein_coding_genes)]
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
    results = process_args()
    seqs = SeqIO.parse(open(results.pb_fasta), 'fasta')
    orfs = pd.read_csv(results.orfs, sep = '\t')
    
    # Filter ORFS based on score and whether protein coding
    orfs = filter_orf_scores(orfs, results.cutoff)
    protein_coding_only = string_to_boolean(results.protein_coding_only)
    if protein_coding_only:
        orfs = filter_protein_coding(orfs, results.protein_coding_file)
    
    pb_seqs,redundant_accs = get_accession_seqs(seqs)
    pb_pseqs = combine_by_sequence(orfs, pb_seqs)
    
    with open(results.red, 'w') as f:
        for acc in redundant_accs:
            f.write(f"{acc}\n")
    
    # Save combined results
    with open(results.combined_tsv, 'w') as ofile, open(results.combined_fasta, 'w') as ofile2:
        ofile.write('protein_sequence\tpb_accs\n')
        for seq, accs in pb_pseqs.items():
            accs_sorted = order_pb_acc_numerically(accs)
            accs_str = '|'.join(accs_sorted)
            ofile.write(seq + '\t' + accs_str + '\n')
            ofile2.write('>' + accs_str + '\n' + seq + '\n')
    
    
    pacbio = pd.read_csv(results.combined_tsv, sep = '\t')
    seqs = SeqIO.parse(open(results.combined_fasta), 'fasta')
    
    pacbio = aggregate_results(pacbio, orfs)
    
    
    with open(results.agg_fasta, "w") as ofile:
        for entry in seqs:
            seq = str(entry.seq)
            pb_acc = entry.id
            pb_row = pacbio[pacbio['pb_accs'] == pb_acc].iloc[0]
            base_acc = pb_row['base_acc']
            gene = orfs[orfs['pb_acc'] == base_acc].iloc[0]['gene']
            ofile.write(f">pb|{base_acc}|fullname GN={gene}\n{seq}\n")
    
    pacbio = pacbio[['pb_accs', 'base_acc', 'FL', 'CPM']]
    pacbio.to_csv(results.agg_tsv, sep = '\t', index = False) 
    
    
if __name__ == "__main__":
    main()