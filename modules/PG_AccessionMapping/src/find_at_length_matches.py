#!/usr/bin/env python3

# find at-length matches (see notes) for query to subject alignments
# 'at-length' match criteria
#  - length of query match equals the full length of the database entry
#  - no gaps
#  - percent matched is >95%

from collections import defaultdict

matches = defaultdict(list)  # query acc -> list of subject acc with match
info = defaultdict(list)  # (isoacc_6k, gene_ensp) -> [pident, gaps]
un_id_w_at_len_match = set()
un_id_all = set()
for line in open('../1_uniprot_to_gencode_blast/run_UN_to_GC_blast_result.tsv').readlines():
    fields = line.split('\t')
    un_id = fields[0]
    un_len = fields[1]
    gc_id = fields[4]
    gc_len = fields[5]
    match_len = fields[10]
    pident = float(fields[11])/float(match_len)
    gaps = int(fields[12])
    # Determine if orf an 'at-length' match.
    if (un_len == gc_len == match_len) and pident > 0.95 and gaps == 0:
        matches[un_id].append((gc_id, match_len, pident))
        un_id_w_at_len_match.add(un_id)
    un_id_all.add(un_id)

# uniprot with no at-length match
un_id_no_match = list(un_id_all - un_id_w_at_len_match)

# write out best alignment for uniprot acc without an at-lenght match with gencode
with open('uniprot_no_at_len_match_to_gencode.tsv', 'w') as ofile:
    ofile.write('uniprot_gene\tuniprot_acc\tuniprot_len\tgencode_gene\tgencode_acc\tgencode_len\tmatch_len\tpercent_identity\tnum_gaps\n')
    for line in open('../1_uniprot_to_gencode_blast/run_UN_to_GC_blast_result.tsv').readlines():
        fields = line.split('\t')
        un_id = fields[0]
        if un_id in un_id_no_match:
            un_id_no_match.remove(un_id)
            un_len = fields[1]
            gc_id = fields[4]
            gc_len = fields[5]
            match_len = fields[10]
            pident = float(fields[11])/float(match_len)
            gaps = int(fields[12])
            uniprot_gene = un_id.split('|')[2].split('_HUMAN')[0]
            gencode_gene = gc_id.split('|')[6]
            ofile.write('\t'.join(map(str,[uniprot_gene, un_id, un_len, gencode_gene, gc_id, gc_len, match_len, pident, gaps])) + '\n')

# write out uniprot with no blast result against gencode
# result - all uniprot accessions blasted to at least one gencode entry
with open('uniprot_no_blast_to_gencode.tsv', 'w') as ofile:
    ofile.write('\n'.join(sorted(un_id_no_match)))


with open('uniprot_at_len_matches_to_gencode.tsv', 'w') as ofile, open('multimapping_cases.tsv', 'w') as ofile2:
    ofile.write('uniprot_gene\tuniprot_acc\tgencode_gene\tgencode_acc\tmatch_len\tpident\n')
    ofile2.write('These are cases where a uniprot sequence mapped well to multiple gencode entries\n')
    for uniprot_acc, gc_info in matches.items():
        uniprot_gene = uniprot_acc.split('|')[2].split('_HUMAN')[0]
        # take the best matching entry, even if there are multiple mappings
        gencode_acc, match_len, pident = gc_info[0]
        if len(gc_info) > 1:
            ofile2.write(str(gc_info) + '\n')
        gencode_gene = gencode_acc.split('|')[6]
        odata = [uniprot_gene, uniprot_acc, gencode_gene, gencode_acc, match_len, pident]
        ofile.write('\t'.join(map(str, odata)) + '\n')
