# %%

## Compile Key Transcriptomic Information ##
# Prepare two tables:
# 1. Sqanti Isoform Table 
# 2. Gene Based Info Tabe 

#### Import Modules ####
# Python Modules #
import numpy as np
import pandas as pd 
import os

#### Input Files ####
sqanti_out = '../../data/jurkat_classification.txt'
tpm_file =  '../../data/jurkat_gene_kallisto.tsv'
ribodep_tpm = '../../data/kallist_table_rdeplete_jurkat.tsv' # expects normalized data
ensg_to_gene = "../../results/PG_ReferenceTables/ensg_to_gene.tsv"
enst_to_isoname = "../../results/PG_ReferenceTables/enst_to_isoname.tsv"
gene_len_stats_tab =  '../../results/PG_ReferenceTables/gene_len_stats.tsv'

#### Outputs ###
# 1. Sqanti_isoform_table
# 2. Gene_level_info_table

#### Part 1 : Prepare Isoform Information Table from sqanti Output ####

def sqtab(sqanti_out, ensg_to_gene, enst_to_isoname):
    """
    Sorts data from Sqanti output 
    """

    # Import Sqanti Output File
    cols = ['isoform', 'length', 'structural_category','associated_gene','associated_transcript','subcategory', 'FL'] 
    data = pd.read_csv(sqanti_out, delimiter="\t", usecols = cols)
    data.columns = ['pb_acc', 'len', 'cat', 'gene','transcript', 'cat2', 'fl_cts']

    # Convert Structural Categories to Acronyms
    data.replace({"novel_not_in_catalog":"NNC","novel_in_catalog":"NIC","incomplete-splice_match":"ISM","full-splice_match":"FSM"}, inplace = True)

    # Filter out any cat that is not FSM, ISM, NNC or NIC
    f= ['FSM', 'ISM', 'NNC', "NIC"]
    fdata = data[data.cat.isin(f)]

    # Normalize fl_cts to cpm 
    sum = fdata['fl_cts'].sum(skipna=True)
    fdata['cpm'] = 1000000*fdata['fl_cts']/sum

    ## Finding and Replacing Gene Information ##
    # Import Human Readable Gene Info and rename columns
    gen_name = pd.read_csv(ensg_to_gene, delimiter="\t", header=None)
    gen_name.columns = ['A', 'B']

    # Make a Dictionary with Columns A and B
    gdict = pd.Series(gen_name.B.values,index=gen_name.A).to_dict()

    # Use Dictionary to Find Gene Names and Replace them with Human Readable Genes
    df = fdata[['gene']]
    fdata['gene'] = fdata['gene'].map(gdict).fillna(fdata['gene'])
    fdata.drop(fdata[fdata['gene'] == df['gene']].index, inplace=True)
  
    ## Finding and Replacing Transcript Information ##
    # Import Human Readable Transcript Info and rename columns
    trans = pd.read_csv(enst_to_isoname, delimiter="\t", header=None)
    trans.columns = ['A', 'B']

    # Make a Dictionary and Replace Transcript Names
    tdict = pd.Series(trans.B.values, index=trans.A).to_dict()
    fdata['transcript'] = fdata['transcript'].map(tdict).fillna(fdata['transcript'])

    return fdata

    print("Isoform Table from sqanti output has been prepared")


def abund(sq_isotab, tpm_file):
    """
    Prepare a gene, cpm, tpm table from sqanti and kallisto output
    """

    # Sort CPM Data
    cpm_data = sq_isotab[['gene', 'cpm']]
    cpm_by_gene = cpm_data.groupby(['gene']).agg(cpm = ('cpm', 'sum')).reset_index(level=['gene'])

    # Sort Kallisto TPM Data
    tpm_by_gene = pd.read_csv(tpm_file, delimiter='\t')
    tpm_by_gene['gene'] = tpm_by_gene['gene'].str.replace('-', '_')

    # Merge
    ab = pd.merge(cpm_by_gene, tpm_by_gene, how='right', on='gene')
    return ab

#-------------------------------------------------------------------------------------------------------
# Main 

# If results folder does not exist, make it
rdir = '../../results/LR_TranscriptomeSummary'
if not os.path.exists(rdir):
    os.mkdir(rdir)

# Make Sqanti Table 
sq_isotab = sqtab(sqanti_out, ensg_to_gene, enst_to_isoname)

# Write Sqanti Dataframe as TSV File
sq_isotab.to_csv("../../results/LR_TranscriptomeSummary/sqanti_isoform_info.tsv", sep="\t", index= False, na_rep='0')

# Make Abundance Table and Merge with Gene_Length_Stats Table 
ab_tab = abund(sq_isotab, tpm_file)
gene_len_stats = pd.read_csv(gene_len_stats_tab, sep = '\t')
gen_lenab = pd.merge(gene_len_stats, ab_tab, how="right", on='gene')

# Make and Merge with PolyA Table 
ribo = pd.read_csv(ribodep_tpm, sep='\t')
rgen = ribo.groupby(['gene']).agg(rtpm=('tpm', 'sum')).reset_index()

# option to output log ratio for the ribosomal data
#rgen['log(rtpm+1)'] = np.log10(rgen['rtpm'] + 1)
#ab_tab['log(tpm+1)'] = np.log10(ab_tab['tpm'] + 1)

ratio = pd.merge(rgen, ab_tab, how = 'left', on='gene')
ratio['rtpm/tpm'] = ratio['rtpm']/ratio['tpm']
ratio_tab = ratio[['gene', 'rtpm/tpm']]
gen_tab = pd.merge(gen_lenab, ratio_tab, how='left', on='gene')

# Output Table 
gen_tab.to_csv("../../results/LR_TranscriptomeSummary/gene_level_tab.tsv", sep="\t", index= False, na_rep='0')

# %%