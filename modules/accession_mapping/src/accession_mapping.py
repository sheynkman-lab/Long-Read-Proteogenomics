#!/usr/bin/env python3

# given protein fasta files for gencode, uniprot, and pacbio databases, map accessions based on sequence similarity
# two sequences are considered a match if both protein sequences match end-to-end (same length) and contain no more than 2 AA mismatches
# however, mapping will first be attempted as perfect matches or 1 AA mismatch
# assumes that the two input fastas contain unique sequence and the representative acc is listed

# %%

from Bio import SeqIO
from collections import defaultdict
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('--gencode_fasta', action='store', dest='gencode_fasta')
parser.add_argument('--uniprot_fasta', action='store', dest='uniprot_fasta')
parser.add_argument('--pacbio_fasta',action='store',dest='pacbio_fasta')
results = parser.parse_args()

# for testing
# results = argparse.Namespace(gencode_fasta='./gencode.fasta', uniprot_fasta='./uniprot_reviewed_canonical_and_isoform.fasta', pacbio_fasta='./jurkat_orf_aggregated.fasta')

# %%

## populate dicitionaries containing exact matches

def num_db_represented(accs):
    num_dbs = len(set([acc.split('|')[0] for acc in accs]))
    return num_dbs

# prot_seq -> <list of gc/un/pb accessions>
seqs = defaultdict(list)
for fasta_file in [results.gencode_fasta, results.uniprot_fasta, results.pacbio_fasta]:
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        seqs[rec.seq].append(rec.id)

multi_match = defaultdict()
single_match = defaultdict()

for seq, accs in seqs.items():
    if num_db_represented(accs) > 1:
        multi_match[seq] = accs
    elif num_db_represented(accs) == 1:
        single_match[seq] = accs[0]
    else:
        print('error:' + seq + ':' + ','.join(accs))

# %%

## for entries without an exact match to other accessions, attempt to find
## "at-length" matches (same length, up to 1 or 2 mismatches)

def is_at_length_match(seq1, seq2, num_mismatches=1):
    """
    Given two sequences as AA strings of the same length,
    determine if they represent an at-length match.

    Parameters
    -----
    seq1 - string, the first AA sequence
    seq2 = string, the second AA sequence
    """
    mismatches = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    if mismatches <= num_mismatches:
        return True
    else:
        return False

def make_len_to_seq_dict(multi_match):
    """ Make a dictionary of sequence length to list of sequences. """
    len_seqs = defaultdict(list)
    for seq, accs in multi_match.items():
        len_seqs[len(seq)].append(seq)
    return len_seqs

len_seqs = make_len_to_seq_dict(multi_match)
still_single_match = defaultdict(list)
for seq, accs in single_match.items():
    found_an_at_length_match = False
    candidate_at_len_seqs = len_seqs[len(seq)]
    for candidate_seq in candidate_at_len_seqs:
        # try to find if there is a match with one allowable mismatches
        # (note - only first at-length match considered)
        if is_at_length_match(seq, candidate_seq, num_mismatches=1):
            # got a hit, add entry to multi_match pool
            multi_match[candidate_seq].extend(accs)
            found_an_at_length_match = True
            break
    if not found_an_at_length_match:
        # try to find if there is a match with two allowable mismatches
        for candidate_seq in candidate_at_len_seqs:
            if is_at_length_match(seq, candidate_seq, num_mismatches=2):
                # got a hit, add entry to multi_match pool
                multi_match[candidate_seq].extend(accs)
                found_an_at_length_match = True
                break
    if not found_an_at_length_match:
        # at this point, found no at-length matches, put in single pool
        still_single_match[seq] = accs

# %%

## write out mapping results

def get_acc_str(list_of_acc, prefix):
    """ Get acc from certain db. Ouput accession is "stripped" down to isoname. """
    sub_accs = []
    for acc in list_of_acc:
        if acc.startswith(prefix):
            acc_stripped = acc.split('|')[1]
            sub_accs.append(acc_stripped)
    return ','.join(sub_accs)

with open('./accession_map_gencode_uniprot_pacbio.tsv', 'w') as ofile:
    ofile.write('gencode_acc\tuniprot_acc\tpacbio_acc\tseq\n')
    for seq, accs in multi_match.items():
        gc = get_acc_str(accs, 'gc')
        un = get_acc_str(accs, 'sp')
        pb = get_acc_str(accs, 'pb')
        ofile.write('\t'.join([gc, un, pb, str(seq)]) + '\n')
    for seq, accs in still_single_match.items():
        gc = get_acc_str([accs], 'gc')
        un = get_acc_str([accs], 'sp')
        pb = get_acc_str([accs], 'pb')
        ofile.write('\t'.join([gc, un, pb, str(seq)]) + '\n')

with open('./accession_map_stats.tsv', 'w') as ofile:
    ofile.write('number of entries with multiple mappings between')
print(len(multi_match))
print(len(single_match))
print(len(still_single_match))

# %%

## write out mapping summary statistics

acc_map = pd.read_table('./accession_map_gencode_uniprot_pacbio.tsv').iloc[:, 0:3]
mapping_frequencies = acc_map.notnull().groupby(['gencode_acc', 'uniprot_acc', 'pacbio_acc']).size().reset_index()
mapping_frequencies.rename(columns={0: 'number_of_entries'}, inplace=True)

mapping_frequencies.to_csv('./accession_map_stats.tsv', sep='\t', index=None)











