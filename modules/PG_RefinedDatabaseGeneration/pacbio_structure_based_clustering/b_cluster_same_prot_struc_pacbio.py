# cluster pacbio proteins with same isoform structure
# cluster orfs predicted from pacbio transcripts, based on genomic structure

# input - gtf file with cds features denoting the ORFs
# output - gtf with same-orf-structure pacbio entries

from collections import defaultdict


# read in the orf coords from gtf file
pb_coords = defaultdict(lambda: ['', '', []]) # pb_acc -> [chr, strand, coords]
i = 0
for line in open('a_gtf_w_cds.gtf'):
    i += 1
    if i > 10000: break
    wds = line.strip().split('\t')
    if wds[2] == 'CDS':
        chr = wds[0]
        start = wds[3]
        end = wds[4]
        strand = wds[6]
        acc = wds[8]
        pb_coords[acc][0] = chr
        pb_coords[acc][1] = strand
        pb_coords[acc][2].append([int(start), int(end)])

# %%

print(pb_coords)

# %%
# group pb acc by same coord
def make_coord_string(infos):
    # return a coordinate string
    # e.g., input [chr1, +, [[50, 100], [150, 200]]]
    #       output chr1+:50-100,150-200
    chr, strand, coords = infos
    coords_sorted = sorted(coords)
    coord_str = chr + strand + ':' + ','.join(str(start) + '-' + str(end) for start, end in coords_sorted)
    return coord_str

def make_internal_coord_string(infos):
    # return a coordinate string that represents the junctions
    # in the case of mono-exon transcripts, return the start and end
    # e.g., input [chr1, +, [[50, 100], [150, 200]]]
    #       output chr1+:100,150
    chr, strand, coords = infos
    coords_sorted = sorted(coords)
    if len(coords_sorted) == 1:
        # for mono-exonic entries, return the full coords (start/end of exon)
        return make_coord_string(infos)
    else:
        all_coords = ','.join(str(start) + '-' + str(end) for start, end in coords_sorted)
        # trim the coords
        internal_coords = '-'.join(all_coords.split('-')[1:-1])
        internal_coords_str = chr + strand + ':' + internal_coords
        return internal_coords_str

# group pb accessions by coordinate structures

coord_accs = defaultdict(lambda: list()) # coord_acc -> list of pb_acc_w_ct
# e.g., chr1+:150-200,250-300 -> PB.1.1_10
icoord_accs = defaultdict(lambda: list()) # internal_coords_acc -> pb_acc_w_ct
# e.g., chr1+:200,250 -> PB.1.1_10

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

for acc, infos in pb_coords.items():
    coord_str = make_coord_string(infos)
    coord_accs[coord_str].append(acc)

    icoord_str = make_internal_coord_string(infos)
    icoord_accs[icoord_str].append(acc)

with open('b_pacbio_same_prot_struc_clusters.tsv', 'w') as ofile:
    ofile.write('protein_coordinates\tpb_accs\n')
    for coord_str, accs in coord_accs.items():
        ordered_accs = order_pb_acc_numerically(accs)
        accs_str = '|'.join(ordered_accs)
        ofile.write(coord_str + '\t' + accs_str + '\n')

with open('b_pacbio_same_internal_junc_prot_struc_clusters.tsv', 'w') as ofile:
    ofile.write('protein_coordinates\tpb_accs\n')
    for coord_str, accs in icoord_accs.items():
        ordered_accs = order_pb_acc_numerically(accs)
        accs_str = '|'.join(ordered_accs)
        ofile.write(coord_str + '\t' + accs_str + '\n')
