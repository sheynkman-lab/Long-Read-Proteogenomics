#!/usr/bin/env python3
#%%
import pandas as pd 
import argparse

def classify_protein_splice_fsm(row):
    if row.pr_splice_cat=='full-splice_match' and row.pr_splice_subcat=='multi-exon':
        ## if nterm and cterm match
        if row.pr_nterm_diff==0 and row.pr_cterm_diff==0:
            return 'pFSM,known_nterm_known_splice_known_cterm'
        ## if nterm and cterm do not match
        elif (row.pr_nterm_diff!=0 or row.pr_cterm_diff!=0) and row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
            return 'pNIC,combo_nterm_cterm'
        ## nterm is matched, cterm is not matched
        elif row.pr_nterm_diff==0 and row.pr_cterm_diff!=0:
            return 'pNNC,known_nterm_known_splice_novel_cterm'
        ## nterm is not matched, cterm is matched
        ## all cases below are when row.pr_nterm_gene_diff!=0
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff==0:
            if row.pr_nhang>0:
                return 'pNNC,novel_nterm_known_splice_known_cterm'
            elif row.pr_nhang<=0:
                if row.utr_cat=='subset':
                    return 'pISM,ntrunc'
                elif row.utr_cat=='unique':
                    return 'pNNC,novel_nterm_known_splice_known_cterm'
        ## nterm is not matched, cterm is not matched
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff!=0:
            if row.pr_nterm_gene_diff==0:
                return 'pNNC,known_nterm_known_splice_novel_cterm'
            elif row.pr_nterm_gene_diff!=0:
                if row.utr_cat=='subset':
                    return 'pISM,ntrunc'
                elif row.utr_cat=='unique':
                    if row.pr_cterm_gene_diff==0:
                        return 'pNNC,novel_nterm_known_splice_known_cterm'
                    elif row.pr_cterm_gene_diff!=0:
                        return 'pNNC,novel_nterm_known_splice_novel_cterm'
        else:
            return 'orphan_fsm'
    return ''

def classify_protein_splice_ism(row):
    if row.pr_splice_cat=='incomplete-splice_match' and 'mono-exon' not in row.pr_splice_subcat:
        ### cases where nterm matches but cterm differs
        if row.pr_nterm_diff==0 and row.pr_cterm_diff!=0:
            if row.pr_cterm_gene_diff!=0:
                return 'pNNC,known_nterm_known_splice_novel_cterm'
            elif row.pr_cterm_gene_diff==0:
                return 'pNIC,combo_nterm_cterm'
        ### cases where nterm and cterm match (retained intron causes unique splicing)
        elif row.pr_nterm_diff==0 and row.pr_cterm_diff==0:
            return 'pNIC,known_nterm_combo_splice_known_cterm'
        ### cases where nterm differs, but cterm matches
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff==0:
            if row.pr_nterm_gene_diff==0:
                return 'pNIC,combo_nterm_cterm'
            else: # row.pr_nterm_gene_diff!=0
                if row.pr_nhang<=0:
                    # determine if nterm is valid novel nterm or an ntrunc
                    if row.utr_cat=='subset':
                        return 'pISM,ntrunc'
                    elif row.utr_cat=='unique':
                        return 'pNNC,novel_nterm_known_splice_known_cterm'
                    else:
                        raise Exception('Invalid utr_cat - needs to be subset or unique')
                elif row.pr_nhang>0:
                    return 'pNNC,novel_nterm_known_splice_known_cterm'
        ### cases where nterm and cterm diff
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff!=0:
            if row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
                return 'pNIC,combo_nterm_cterm'
            elif row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff!=0:
                return 'pNNC,known_nterm_known_splice_novel_cterm'
            elif row.pr_nterm_gene_diff!=0:
                if row.utr_cat=='subset':
                    return 'pISM,ntrunc'
                elif row.utr_cat=='unique':
                    if row.pr_cterm_gene_diff==0:
                        return 'pNNC,novel_nterm_known_splice_known_cterm'
                    elif row.pr_cterm_gene_diff!=0:
                        return 'pNNC,novel_nterm_known_splice_novel_cterm'
        else:
            return 'orphan_ism'
    return ''

def classify_protein_splice_nic(row):
    if row.pr_splice_cat=='novel_in_catalog' and 'mono-exon' not in row.pr_splice_subcat:
        if row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
            return 'pNIC,known_nterm_combo_splice_known_cterm'
        elif row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff!=0:
            return 'pNNC,known_nterm_combo_splice_novel_cterm'
        elif row.pr_nterm_gene_diff!=0:
            if row.utr_cat=='subset':
                return 'pISM,ntrunc'
            elif row.utr_cat=='unique':
                if row.pr_cterm_gene_diff==0:
                    return 'pNNC,novel_nterm_combo_splice_known_cterm'
                else:
                    return 'pNNC,novel_nterm_combo_splice_novel_cterm'
        else:
            return 'orphan_nic'
    return ''

def classify_protein_splice_nnc(row):
    if row.pr_splice_cat=='novel_not_in_catalog':
        if row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
            return 'pNNC,known_nterm_novel_splice_known_cterm'
        elif row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff!=0:
            return 'pNNC,known_nterm_novel_splice_novel_cterm'
        elif row.pr_nterm_gene_diff!=0:
            if row.utr_cat=='subset':
                return 'pISM,ntrunc' 
            elif row.utr_cat=='unique':
                if row.pr_cterm_gene_diff==0:
                    return 'pNNC,novel_nterm_novel_splice_known_cterm'
                else:
                    return 'pNNC,novel_nterm_novel_splice_novel_cterm'
        else:
            return 'orphan_nnc'
    return ''

def classify_protein_splice_monoexon(row):
    if row.pr_splice_subcat=='mono-exon' or row.pr_splice_subcat=='mono-exon_by_intron_retention':
        if row.pr_splice_cat=='full-splice_match' and row.pr_nterm_diff==0 and row.pr_cterm_diff==0:
            return 'pFSM,mono-exon'
        elif row.pr_splice_cat=='intergenic':
            return 'intergenic,mono-exon'
        # ignoring previous trunc/altnterm
        else:
            return 'orphan_monoexon,mono-exon'
    return ''

def classify_protein_splice_misc(row):
    if row.pr_splice_subcat=='multi-exon':
        if row.pr_splice_cat=='intergenic':
            return 'intergenic,multi-exon'
        elif row.pr_splice_cat=='genic':
            return 'genic,multi-exon'
        elif row.pr_splice_cat=='antisense':
            return 'antisense,multi-exon'
        elif row.pr_splice_cat=='fusion':
            return 'fusion,multi-exon'
    return ''


def classify_protein(row):
    fsm_classiciation = classify_protein_splice_fsm(row)
    ism_classification = classify_protein_splice_ism(row)
    nic_classification = classify_protein_splice_nic(row)
    nnc_classification = classify_protein_splice_nnc(row)
    mono_classification = classify_protein_splice_monoexon(row)
    misc_classification = classify_protein_splice_misc(row)
    agg_classification = fsm_classiciation + ism_classification + nic_classification + nnc_classification + mono_classification + misc_classification 
    return (
        agg_classification
    )
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sqanti_protein', action='store', dest='sqanti_protein', help='sqanti protein classification file')
    parser.add_argument('--name', action='store',dest='name')
    parser.add_argument('--dest_dir', action='store',dest='dest_dir',default='./')
    args = parser.parse_args()

    protein_classification = pd.read_table(args.sqanti_protein)
    protein_classification['protein_classification'] = protein_classification.apply(classify_protein, axis = 1)
    cats = protein_classification["protein_classification"].str.split(",", n = 1, expand = True)
    protein_classification['protein_classification_base'] = cats[0]
    protein_classification['protein_classification_subset'] = cats[1]

    protein_classification.to_csv(f'{args.dest_dir}{args.name}.protein_classification.tsv', sep='\t', index = False)

#%%
if __name__ == "__main__":
    main()
# %%
