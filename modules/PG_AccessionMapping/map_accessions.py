# given two protein fasta files, map the accessions based on matching sequences
# two sequences are considered a match if both protein sequences match end-to-end (same length) and contain no more than 2 AA mismatches
# assumes that the two input fastas contain unique sequence and the representative acc is listed

# %%

from Bio import SeqIO
from collections import defaultdict


gc_fpath = '../../result/gencode.fasta'
pb_fpath = '../../result/protein_coding_only/jurkat.collapsed_orf_aggregated.fasta'



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


gc = read_fasta_into_dict(gc_fpath)
pb = read_fasta_into_dict(pb_fpath)

len_pbs = make_len_to_seq_list_dict(pb)

matches = defaultdict(list) # gc isoname -> list of pb matches
for isoname, gcseq in gc.items():
    pb_len_matches = len_pbs[len(gcseq)]
    for (acc, pbseq) in pb_len_matches:
        if is_at_length_match(gcseq, pbseq):
            matches[isoname].append(acc)

matches_perfect = defaultdict(list) # gc isoname -> list of pb matches
for isoname, gcseq in gc.items():
    pb_len_matches = len_pbs[len(gcseq)]
    for (acc, pbseq) in pb_len_matches:
        if is_at_length_match(gcseq, pbseq):
            matches_perfect[isoname].append(acc)

list_of_pbs = list(pb.keys())
list_of_pbs_perfect = list(pb.keys())

with open('gc_pb_map_atlenseq.tsv', 'w') as ofile:
    ofile.write('gencode\tpacbio_atlen\tpacbio_perfect\n')
    for isoname, gcseq in gc.items():
        pb_matches = matches[isoname]
        pb_matches_perfect = matches_perfect[isoname]
        ofile.write(isoname + '\t' + ','.join(pb_matches) + '\t' + ','.join(pb_matches_perfect) + '\n')
        for pb_match in pb_matches:
            if pb_match in list_of_pbs:
                list_of_pbs.remove(pb_match)
    for acc in list_of_pbs:
        ofile.write(f'\t{acc}\n')
    






    


# %%



        
