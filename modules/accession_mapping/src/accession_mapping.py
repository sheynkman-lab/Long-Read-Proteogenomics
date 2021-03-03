# given two protein fasta files, map the accessions based on matching sequences
# two sequences are considered a match if both protein sequences match end-to-end (same length) and contain no more than 2 AA mismatches
# assumes that the two input fastas contain unique sequence and the representative acc is listed

# %%

from Bio import SeqIO
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref_fasta', action='store', dest='ref_fasta')
parser.add_argument('--other_fasta',action='store',dest='other_fasta')
parser.add_argument('--ref_name', action='store',dest='ref_name')
parser.add_argument('--other_name', action='store',dest='other_name')

# gc_fpath = '../../result/gencode.fasta'
# pb_fpath = '../../result/protein_coding_only/jurkat.collapsed_orf_aggregated.fasta'



def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def is_same_length(seq1, seq2):
    if len(seq1) == len(seq2):
        return True
    return False

def is_at_length_match(seq1, seq2):
    """
    Given two sequences as AA strings, determine if they represent an at-length match.

    Parameters
    -----
    seq1 - string, the first AA sequence
    seq2 = string, the second AA sequence
    """
    mismatches = hamming_distance(seq1, seq2)
    if mismatches <= 2:
        return True
    else:
        return False

def read_fasta_into_dict(fpath):
    seq_dict = defaultdict()
    for rec in SeqIO.parse(fpath, 'fasta'):
        acc = rec.id.split('|')[1]
        seq_dict[acc] = str(rec.seq)
    return(seq_dict)

def make_len_to_seq_list_dict(pb):
    """ Make a dictionary of sequence lenght to entries. """
    len_pb = defaultdict(list) # length -> [(pbacc, seq), (pbacc, seq)]
    for acc, pbseq in pb.items():
        pblen = len(pbseq)
        len_pb[pblen].append((acc, pbseq))
    return(len_pb)


ref_seq_dict = read_fasta_into_dict(results.fasta1)
other_seq_dict = read_fasta_into_dict(results.fasta2)

other_len_dict = make_len_to_seq_list_dict(other_seq_dict)

matches = defaultdict(list) # gc isoname -> list of pb matches
for ref_id, seq in ref_seq_dict.items():
    other_len_matches = other_len_dict[len(seq)]
    for (acc, other_seq) in other_len_matches:
        if is_at_length_match(seq, other_seq):
            matches[ref_id].append(acc)

matches_perfect = defaultdict(list) # gc isoname -> list of pb matches
for ref_id, seq in ref_seq_dict.items():
    other_len_matches = other_len_dict[len(seq)]
    for (other_id, other_seq) in other_len_matches:
        if is_at_length_match(seq, other_seq):
            matches_perfect[ref_id].append(other_id)

list_of_others = list(other_seq_dict.keys())
list_of_others_perfect = list(other_seq_dict.keys())

with open(f'{results.ref_name}_{results.other_name}_map_atlenseq.tsv', 'w') as ofile:
    ofile.write(f'{results.ref_name}\t{results.other_name}_atlen\t{results.other_name}_perfect\n')
    for ref_id, seq in ref_seq_dict.items():
        other_matches = matches[ref_id]
        other_matches_perfect = matches_perfect[ref_id]
        ofile.write(ref_id + '\t' + ','.join(other_matches) + '\t' + ','.join(other_matches_perfect) + '\n')
        for other_match in other_matches:
            if other_match in list_of_others:
                list_of_others.remove(other_match)
    for acc in list_of_others:
        ofile.write(f'\t{acc}\n')
    






    


# %%



        
