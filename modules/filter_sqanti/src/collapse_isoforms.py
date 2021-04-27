#!/usr/bin/env python3
# collapse entries, based on cdna cupcake
#%%
from cupcake.io import GFF
from collections import defaultdict
from cupcake.tofu import compare_junctions
import os
import argparse
from Bio import SeqIO
### helper functions for main (below)

def can_merge(m, r1, r2, internal_fuzzy_max_dist):
    if m == 'subset':
        r1, r2 = r2, r1 #  rotate so r1 is always the longer one
    if m == 'super' or m == 'subset':
        n2 = len(r2.ref_exons)
        if r1.strand == '+':
            # if r2 is monoexonic, it can start after the last exon of r1's last exon
            # if r2 is multiexonic, the last start must be pretty close (fuzzy allowed)
            if n2==1: # r2 is mono-exonic
                return r1.ref_exons[-1].start - r2.ref_exons[-1].start <= internal_fuzzy_max_dist 
            else: return abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[-n2].start <= r2.ref_exons[0].start < r1.ref_exons[-n2].end
        else:
            if n2==1: return r1.ref_exons[0].end - r2.ref_exons[0].end >= -internal_fuzzy_max_dist
            else: return abs(r1.ref_exons[0].end - r2.ref_exons[0].end) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[n2-1].start <= r2.ref_exons[-1].end < r1.ref_exons[n2].end

def filter_out_subsets(recs, internal_fuzzy_max_dist):
    # recs must be sorted by start becuz that's the order they are written
    i = 0
    while i < len(recs)-1:
        no_change = True
        j = i + 1
        while j < len(recs):
            if recs[j].start > recs[i].end: 
                break
            recs[i].segments = recs[i].ref_exons
            recs[j].segments = recs[j].ref_exons
            m = compare_junctions.compare_junctions(recs[i], recs[j], internal_fuzzy_max_dist)
            if can_merge(m, recs[i], recs[j], internal_fuzzy_max_dist):
                if m == 'super': # pop recs[j] 
                    recs.pop(j)
                else:
                    recs.pop(i)
                    no_change = False
            else:
                j += 1
        if no_change: i += 1

def modify_gff_file(sqanti_gtf, output_filename):
    # make a version of the gff where gene_id is pb accession, so compatible with gff reader
# this file 'jurkat_with_cds_clean_acc_for_ucsc_all_chr' does not get read into gff reader correctly
# is there an issue with the exon ranges?
# for now, process the output from sqanti, corrected gtf
    if not os.path.exists(output_filename):
        with open(output_filename, 'w') as ofile:
            for line in open(sqanti_gtf):
                pb_acc = line.split('transcript_id "')[1].split('"')[0]
                pb_locus_acc = 'PB.' + pb_acc.split('.')[1]
                wds = line.split('\t')
                if wds[2] in ('transcript', 'exon'):
                    prefix = wds[0:8]
                    acc_line = 'gene_id "{}"; transcript_id "{}";'.format(pb_acc, pb_locus_acc)
                    ofile.write('\t'.join(prefix + [acc_line]) + '\n')

#%%
def collapse_isoforms(gff_filename, fasta_filename, name):
    
    recs = defaultdict(lambda: [])
    reader = GFF.collapseGFFReader(gff_filename)
    for record in reader:
        assert record.seqid.startswith('PB.')
        pb_cluster = f'PB.{record.seqid.split(".")[1]}'
        recs[pb_cluster].append(record)
    # collapse and write gff
    good = []
    output_gtf_filename = f'{name}_corrected.5degfilter.gff'
    output_gtf_file = open(output_gtf_filename, 'w')
    for pb_cluster, isoforms in recs.items():
        fuzzy_junc_max_dist = 0
        filter_out_subsets(isoforms, fuzzy_junc_max_dist)
        for record in isoforms:
            GFF.write_collapseGFF_format(output_gtf_file, record)
            good.append(record.seqid)
    output_gtf_file.close()
    # write filtered fata file
    isoform_seqs = []
    for record in SeqIO.parse(fasta_filename, 'fasta'):
        if record.id in good:
            isoform_seqs.append(record)
    SeqIO.write(isoform_seqs, f'{name}_corrected.5degfilter.fasta','fasta')
    
#%%
###################################
### main - find 5' deg products ###
###################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', '-oc',action='store', dest= 'name',)
    parser.add_argument('--sqanti_gtf',action='store',dest='sqanti_gtf')
    parser.add_argument('--sqanti_fasta',action='store',dest='sqanti_fasta')
    args = parser.parse_args()
    collapse_isoforms(args.sqanti_gtf,args.sqanti_fasta, args.name)

if __name__ == "__main__":
    main()


# %%
