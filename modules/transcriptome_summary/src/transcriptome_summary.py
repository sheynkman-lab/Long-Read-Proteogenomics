#!/usr/bin/env python3
"""
This module prepares an isoform classification table and a gene info table 

Inputs:
--------------------------------------------------------------------------
 1. sqanti classification file 
 2. kalliso tpm file 
 3. kalliso ribodepletion tpm file 
 4. ensg_to_gene mapping file 
 5. enst_to_isoname mapping file 
 6. gene_len_stats table 
--------------------------------------------------------------------------

Outputs:
--------------------------------------------------------------------------
1. sqanti isoform table 
2. gene level info table 
--------------------------------------------------------------------------

"""


# Import Modules 
import numpy as np 
import pandas as pd
import argparse
import os 

# Define Functions
def sqtab(sqanti_out, ensg_to_gene, enst_to_isoname):
    """
    Sorts data from Sqanti output 

    Note: We are not considering genes from sqanti output like ENSG00000242861.1_ENSG00000196187.12
    """

    # Import Data
    cols = ['isoform', 'length', 'structural_category','associated_gene','associated_transcript','subcategory', 'FL'] 
    data = pd.read_csv(sqanti_out, delimiter="\t", usecols = cols)
    data.columns = ['pb_acc', 'len', 'cat', 'gene','transcript', 'cat2', 'fl_cts']

    # Map categories to acronyms and filter out anything that is not FSM, ISM, NNC or NIC
    data.replace({"novel_not_in_catalog":"NNC","novel_in_catalog":"NIC","incomplete-splice_match":"ISM","full-splice_match":"FSM"}, inplace = True)
    fdata = data[data.cat.isin(['FSM', 'ISM', 'NNC', "NIC"])]

    # Normalize fl_cts to cpm
    fdata['cpm'] = 1000000*fdata['fl_cts']/fdata['fl_cts'].sum(skipna=True)

    # Map gene -> gene_name
    gen_name = pd.read_csv(ensg_to_gene, delimiter="\t", header=None)
    gdict = pd.Series(gen_name.loc[:,1].values,index=gen_name.loc[:,0]).to_dict()
    df = fdata[['gene']]
    fdata['gene'] = fdata['gene'].map(gdict).fillna(df['gene'])

    # Drop cases like ENSG00000242861.1_ENSG00000196187.12 
    fdata.drop(fdata[fdata['gene'] == df['gene']].index, inplace=True)

    # Map enst -> isoname
    trans = pd.read_csv(enst_to_isoname, delimiter="\t", header=None)
    tdict = pd.Series(trans.loc[:,1].values, index=trans.loc[:,0]).to_dict()
    df2 = fdata[['transcript']]
    fdata['transcript'] = fdata['transcript'].map(tdict).fillna(df2['transcript'])
    print("Isoform Table from sqanti output has been prepared")
    return fdata
    

def abund(sq_isotab, tpm_file):
    """
    Prepare a gene, cpm, tpm table from sqanti and kallisto output
    """ 
    # Sort CPM Data
    cpm_data = sq_isotab[['gene', 'cpm']]
    cpm_by_gene = cpm_data.groupby(['gene']).agg(cpm = ('cpm', 'sum')).reset_index(level=['gene']) 

    # Sort Kallisto TPM Data
    tpm_by_gene = pd.read_csv(tpm_file, delimiter='\t')
    # tpm_by_gene['gene'] = tpm_by_gene['gene'].str.replace('-', '_')

    # Merge
    ab = pd.merge(cpm_by_gene, tpm_by_gene, how='right', on='gene')
    return ab

def main():

    # Main Code 
    parser = argparse.ArgumentParser(description='Process transcriptome related input file locations')
    parser.add_argument('--sq_out', '-s', action='store', dest='sqanti_out', help = 'input : Sqanti Classification output location')
    parser.add_argument('--tpm', '-t', action='store', dest='tpm_file',help='Kallisto TPM file location')
    parser.add_argument('--ribo', '-r', action='store', dest='ribodep_tpm', help='Normalized Kallisto Ribodepletion TPM file location')
    parser.add_argument('--ensg_to_gene', '-gmap', action='store', dest='ensg_to_gene', help='ENSG -> Gene Map file location')
    parser.add_argument('--enst_to_isoname', '-imap', action='store', dest='enst_to_isoname', help='ENST -> Isoname Map file location')
    parser.add_argument('--len_stats', '-l', action='store', dest='gene_len_stats_tab', help='Gene Length Statistics table location')
    parser.add_argument('--odir', '-o', action='store', dest='odir', help='Output Directory', default=None)
    results = parser.parse_args()

    # If results folder does not exist, make it
    odir = results.odir
    if odir is not None and not os.path.exists(odir):
        os.mkdir(odir)
    else:
        odir = ''

    # Make Sqanti Isoform Table and output to a TSV
    sq_isotab = sqtab(results.sqanti_out, results.ensg_to_gene, results.enst_to_isoname)
    # sq_isotab['gene'] = sq_isotab['gene'].str.replace('_','-')
    sq_isotab.to_csv(os.path.join(odir, 'sqanti_isoform_info.tsv'), sep="\t", index= False, na_rep='0')

    # Make PB-Gene reference table
    pb_gene = sq_isotab[['pb_acc','gene']]
    # pb_gene.columns = ['isoform','gene']
    pb_gene = pb_gene.drop_duplicates()
    pb_gene.to_csv(os.path.join(odir, 'pb_gene.tsv'), sep="\t", index= False, na_rep='0')

    # Make Abundance Table and Merge with Gene_Length_Stats Table 
    ab_tab = abund(sq_isotab, results.tpm_file)
    gene_len_stats = pd.read_csv(results.gene_len_stats_tab, sep = '\t')
    gen_lenab = pd.merge(gene_len_stats, ab_tab, how="right", on='gene')

    # Make and Merge with PolyA Table 
    ribo = pd.read_csv(results.ribodep_tpm, sep='\t')
    rgen = ribo.groupby(['gene']).agg(rtpm=('tpm', 'sum')).reset_index()

    # option to output log ratio for the ribosomal data
    #rgen['log(rtpm+1)'] = np.log10(rgen['rtpm'] + 1)
    #ab_tab['log(tpm+1)'] = np.log10(ab_tab['tpm'] + 1)

    ratio = pd.merge(rgen, ab_tab, how = 'left', on='gene')
    ratio['rtpm/tpm'] = ratio['rtpm']/ratio['tpm']
    ratio_tab = ratio[['gene', 'rtpm/tpm']]
    gen_tab = pd.merge(gen_lenab, ratio_tab, how='left', on='gene')

    # Output Table 
    gen_tab.to_csv(os.path.join(odir, 'gene_level_tab.tsv'), sep="\t", index= False, na_rep='0')

if __name__ == "__main__":
    main()