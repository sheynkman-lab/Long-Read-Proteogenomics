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
    parser.add_argument('-is', '--seq', action = 'store', dest='seq', help='PacBio fasta sequence input file location')
    parser.add_argument('-or', '--redundent', action='store', dest='red', help = 'Output redundant accession file location')
    parser.add_argument('-oct', '--ctsv', action='store', dest='combined_tsv', help = 'Output combined tsv file location')
    parser.add_argument('-ocf', '--cfasta', action='store', dest='combined_fasta', help = 'Output combined fasta file location')
    parser.add_argument('-ot', '--aggtsv', action='store', dest='agg_tsv', help = 'Output aggregated tsv file location')
    parser.add_argument('-of', '--aggfasta', action='store', dest='agg_fasta', help = 'Output aggregated fasta file location')
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
    
    fl_dict = pd.Series(orfs.FL.values,index=orf.pb_acc).to_dict()
    cpm_dict = pd.Series(orfs.CPM.values,index=orf.pb_acc).to_dict()
    
    pacbio['FL'] = pacbio['accessions'].apply(lambda accs : get_total(accs, fl_dict))
    pacbio['CPM'] = pacbio['accessions'].apply(lambda accs : get_total(accs, cpm_dict))
    return pacbio
    
    
def main():
    results = process_args()
    seqs = SeqIO.parse(open(results.seq), 'fasta')
    orfs = pd.read_csv(results.orfs, sep = '\t')               
    pb_seqs,redundant_accs = get_accession_seqs(seqs)
    pb_pseqs = aggregate_by_sequence(orfs, pb_seqs)
    
    
    with open(results.red, 'w') as f:
    for acc in redundant_accs:
        f.write(f"{acc}\n")
        
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