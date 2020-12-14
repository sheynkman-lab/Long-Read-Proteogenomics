#!/usr/bin/env python3

# categorization of protein group comparisons across databases


import pandas as pd
import numpy as np
from collections import defaultdict


isomap_toy = pd.read_table('AccessionKey_formated.tsv')
# isomap = pd.read_table('isoform_accession_mapping_table.tsv')

isomap_toy['un'] = isomap_toy['uniprot_acc'].str.split('|').str[1]
isomap_toy['pb'] = isomap_toy['pacbio_acc']
isomap_toy['gc'] = isomap_toy['gencode_acc_y'].str.split('|').str[0]
imap = isomap_toy[['un', 'pb', 'gc']]
imap['idx'] = np.arange(len(imap))


# create pb to ref accession map
pb_ref = defaultdict() # pb_acc -> [pb, un, gc] (no entry if no accession)

def return_list_with_nans_filtered_out(input_list):
    filtered_list = []
    for item in input_list:
        if pd.isnull(item): continue
        filtered_list.append(item)
    return filtered_list

for i, row in imap.iterrows():
    un, pb, gc, idx = row
    acc_list = return_list_with_nans_filtered_out([pb, un, gc])
    pb_ref[pb] = acc_list[-1]


# create un to gc accession map


# %%

gc = pd.read_table('./MockGENCODE_ProteinGroups_formated.tsv')
pb = pd.read_table('./MockPacBio_ProteinGroups_formated.tsv')

gc_acc = gc['Protein Accession'].dropna().to_list()
# are all protein groups input unique?
gc_acc
len(gc_acc) == len(set(gc_acc))

# are all protein groups input unique?
pb_acc = pb['Protein Accession'].dropna().to_list()
pb_acc
len(pb_acc) == len(set(pb_acc))

grps1 = gc_acc
grps2 = pb_acc
unassigned_grps = [grps1, grps2]

# determine which group of protein accessions is pacbio/gencode

def grp_is_derived_from_pacbio(grp):
    peek_acc = grp[0].split('|')[0]
    if peek_acc.startswith('PB.'):
        return True
    else:
        return False

def grp_is_derived_from_gencode(grp):
    peek1 = grp[0].split('|')[0]
    peek2 = grp[1].split('|')[0]
    truth_table = []
    for peek in [peek1, peek2]:
        if '-' not in peek:
            truth_table.append(False)
        else:
            iso_num = int(peek.strip().split('-')[-1])
            if iso_num > 100:
                truth_table.append(True)
            else:
                truth_table.append(False)
    if sum(truth_table) == 2:
        return True
    else:
        return False

def determine_source_database_for_grps(unassigned_grps):
    """From a list of list of protein groups, determine the source database"""
    grps = [] # e.g., ['un': grp, 'pb': grp]
    for grp in unassigned_grps:
        if grp_is_derived_from_pacbio(grp):
            grps.append(['pb':grp])
        elif grp_is_derived_from_gencode(grp):
            grps.append(['gc':grp])
        else
            grps.append(['un':grp])

grps = determine_source_database_for_grps(unassigned_grps)






pb_converted




# %%
