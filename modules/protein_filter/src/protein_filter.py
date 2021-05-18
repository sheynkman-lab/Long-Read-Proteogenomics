#!/usr/bin/env python3
# filter out intergenic, candidate truncates
# start with file ben sent me with the protein classifications
# make a pclass table (clean) for ucsc displayed info

#%%


import pandas as pd
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser()
parser.add_argument('--protein_classification',action='store',dest='protein_classification')
parser.add_argument('--gencode_gtf',action='store',dest='gencode_gtf')
parser.add_argument('--protein_fasta',action='store',dest='protein_fasta')
parser.add_argument('--sample_cds_gtf',action='store',dest='sample_cds')
parser.add_argument('--min_junctions_after_stop_codon',action='store', dest='min_junctions_after_stop_codon',type=int,default=2)
parser.add_argument('--name',action='store',dest='name')

args = parser.parse_args()

#%%


def determine_if_pb_should_be_filtered(row, min_junc_after_stop_codon):
    """PB should be filtered if NMD, a truncation, or protein classification
    is not likely protein coding (intergenic, antisense, fusion,...)

    Args:
        row (pandas Series): protein classification row
        min_junc_after_stop_codon (int): mininum number of junctions after stop
        codon a protein can have. used in NMD determination

    Returns:
        int: 1 if should be filtered, 0 if should not be filtered
    """
    # filter out pbs that are artifacts or noncoding
    pclass = str(row['protein_classification'])
    num_junc_after_stop_codon = int(row['num_junc_after_stop_codon'])
    pclass_base_to_keep = ['pFSM','pNIC']
    pclass_base = str(row['protein_classification_base'])
    if pclass_base not in pclass_base_to_keep and num_junc_after_stop_codon > min_junc_after_stop_codon:
        return 1
    elif 'trunc' in pclass: 
        return 1
    elif 'intergenic' in pclass:
        return 1
    elif 'antisense' in pclass:
        return 1
    elif 'fusion' in pclass:
        return 1
    elif 'orphan' in pclass:
        return 1
    elif 'genic' in pclass:
        return 1
    return 0

def get_short_psubclass_descriptor(psubclass):
    """Generates a short descriptor of the subclass

    Args:
        psubclass (string): protein classification subclass

    Raises:
        Exception: unknown subclass given

    Returns:
        string: short subclass descriptor
    """
    # derive a shorter psubclass for viewing on ucsc browser
    if psubclass == 'known_nterm_novel_splice_known_cterm':
        return 'kn_ns_kc'
    elif psubclass == 'known_nterm_known_splice_known_cterm':
        return 'kn_ks_kc'
    elif psubclass == 'known_nterm_combo_splice_known_cterm':
        return 'kn_cs_kc'
    elif psubclass == 'known_nterm_novel_splice_novel_cterm':
        return 'kn_ns_nc'
    elif psubclass == 'known_nterm_combo_splice_novel_cterm':
        return 'kn_cs_nc'
    elif psubclass == 'known_nterm_known_splice_novel_cterm':
        return 'kn_ks_nc'
    elif psubclass == 'novel_nterm_known_splice_known_cterm':
        return 'nn_ks_kc'
    elif psubclass == 'mono-exon':
        return 'mono'
    elif psubclass == 'ntrunc':
        return 'ntrunc'
    elif psubclass == 'combo_nterm_cterm':
        return 'cnc'
    elif psubclass == 'novel_nterm_known_splice_novel_cterm':
        return 'nn_ks_nc'
    elif psubclass == 'novel_nterm_novel_splice_known_cterm':
        return 'nn_ns_kc'
    elif psubclass == 'novel_nterm_combo_splice_novel_cterm':
        return 'nn_cs_nc'
    elif psubclass == 'novel_nterm_combo_splice_known_cterm':
        return 'nn_cs_kc'
    elif psubclass == 'novel_nterm_novel_splice_novel_cterm':
        return 'nn_ns_nc'
    elif psubclass == 'multi-exon':
        return 'multi'
    elif psubclass == 'ctrunc':
        return 'ctrunc'
    else:
        raise Exception('Invalid psubclass:' + psubclass)

prot = pd.read_table(args.protein_classification)
prot = prot.dropna(subset=['protein_classification'])
prot['filter_status'] = prot.apply(lambda row: determine_if_pb_should_be_filtered(row, args.min_junctions_after_stop_codon), axis=1)


prot['pclass'] = prot['protein_classification'].str.split(',').str[0]
prot['pclass'] = prot['pclass'].fillna('-')
prot['psubclass'] = prot['protein_classification'].str.split(',').str[1]
prot['psubclass_short'] = prot['psubclass'].apply(get_short_psubclass_descriptor)

# output info needed for ucsc track visualization
prot['orf_calling_confidence']
prot['orf_conf'] = prot['orf_calling_confidence'].str[0:3]
# add genename column
# ensgs = {}
# for line in open(args.gencode_gtf):
#     if line.startswith('#'): continue
#     seqname,source,feature,start,end,score,strand,phase,attributes=line.split('\t')
#     if feature.strip()=='gene':
#         ensg = line.split('gene_id "')[1].split('"')[0]
#         gene = line.split('gene_name "')[1].split('"')[0]
#         ensgs[ensg] = gene
# prot['gene'] = prot['pr_gene'].map(ensgs)
# end add genename columnd
(
    prot
        .filter(['pb', 'tx_gene','pr_gene', 'is_nmd', 'has_stop_codon', 'orf_conf', 'CPM', 'pclass', 'psubclass', 'psubclass_short', 'filter_status'])
        .to_csv(f'./{args.name}_w_class_info.tsv', sep='\t', index=None)
)

(
    prot
        .query('filter_status==0')
        .filter(['pb', 'tx_gene','pr_gene', 'is_nmd', 'has_stop_codon', 'orf_conf', 'CPM', 'pclass', 'psubclass', 'psubclass_short',])
        .to_csv(f'./{args.name}.classification_filtered.tsv', sep='\t', index=None)
)
#%%
filtered_proteins = prot.query('filter_status==0')
filtered_accs = set(filtered_proteins['pb'])
filtered_fasta = []
for record in SeqIO.parse(args.protein_fasta, 'fasta'):
    acc = record.description.split('|')[1]
    if acc in filtered_accs:
        filtered_fasta.append(record)

SeqIO.write(filtered_fasta, f'{args.name}.filtered_protein.fasta','fasta')

pb_pclass_dict = pd.Series(filtered_proteins.protein_classification_base.values,index=filtered_proteins.pb).to_dict()

# write filtered cds.gtf file
with open(args.sample_cds, 'r') as ifile, open(f'{args.name}_with_cds_filtered.gtf','w') as cds_ofile:
    for line in ifile:
        seqname,source,feature,start,end,score,strand,phase,attributes=line.split('\t')
        gene_info,pb_acc,cpm = attributes.split("|")
        if pb_acc in filtered_accs:
            # cds_line='\t'.join([seqname,source,feature,start,end,score,strand,phase, attributes])
            # pclass_attributes='|'.join([gene_info,pb_acc,pb_pclass_dict[pb_acc]])
            # pclass_line='\t'.join([seqname,source,feature,start,end,score,strand,phase, pclass_attributes]) + '\n'
            both_attributes='|'.join([gene_info,pb_acc,pb_pclass_dict[pb_acc],cpm])
            both_line=pclass_line='\t'.join([seqname,source,feature,start,end,score,strand,phase, both_attributes])
            cds_ofile.write(both_line)




# %%
