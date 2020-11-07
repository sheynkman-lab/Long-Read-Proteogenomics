# temp scripts for different parts of the pipline

import pandas as pd
from Bio import SeqIO
from Bio import Seq
from collections import defaultdict



# extract protein sequence from pb transcript and given cpat orf ranges

# read in the transcript sequences
pb_seqs = defaultdict() # pb_acc -> transcript_seq
seqs = SeqIO.parse(open('../../a_SQANTI3_out/jurkat_corrected.fasta'), 'fasta')
with open('a_redundant_pb_entries.txt', 'w') as ofile:
    for entry in seqs:
        seq = str(entry.seq)
        pb_acc = entry.id
        # skip duplicates
        if pb_acc in pb_seqs:
            ofile.write(pb_acc + '\n')
            continue
        pb_seqs[pb_acc] = seq

# read in cpat-predicted orf ranges
df = pd.read_table('../0_cpat_analysis_output_all_orfs/jurkat_cpat.ORF_prob.best.tsv')
df = df[['seq_ID', 'ORF_ID', 'ORF_frame', 'ORF_start', 'ORF_end']]
df.columns = ['pb_acc', 'orf_acc', 'start', 'end', 'len']

# extract, translate, and cluster protein sequences
pb_pseqs = defaultdict(lambda: list()) # protein_seq -> list of pb acc
for index, row in df.iterrows():
    pb_acc, _, start, end, olen = list(row)
    seq = pb_seqs[pb_acc]
    orf_seq = seq[start-1:end]
    prot_seq = Seq.translate(orf_seq, to_stop=True)
    pb_pseqs[prot_seq].append(pb_acc)

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


# write out protein cluster data
with open('a_pacbio_protein_clusters.tsv', 'w') as ofile, open('a_pacbio_protein_clusters.fasta', 'w') as ofile2:
    ofile.write('protein_sequence\tpb_accs\n')
    for seq, accs in pb_pseqs.items():
        accs_sorted = order_pb_acc_numerically(accs)
        accs_str = '|'.join(accs_sorted)
        ofile.write(seq + '\t' + accs_str + '\n')
        ofile2.write('>' + accs_str + '\n' + seq + '\n')
