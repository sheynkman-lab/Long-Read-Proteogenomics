# module to generate a table of PacBio isoforms which contains information of relevance to the 
# proteogenomics framework and proteomics interpretation from sqanti outputs.
#
#### READ ME ####
# This function generates a table of pacbio isoforms from sqanti output data
#
# 
# 
# sqtab(sqanti_out, ensg_to_gene, enst_to_trans, tabpath):
#                                               inputs:
#                                                   sqanti_output file path
#                                                   ensg_to_gene file path 
#                                                   enst_to_trans file path 
#                                                   tab_path: path to output table file 
#                                                       default is in current folder
#                                               output: 
#                                                   isoform table file 
#                                                   pandas isoform table 
#
# TO FIX: Some of the sqanti pb_acc map to more than 1 gene. I have removed those genes here, but we need to think about how to fix this


# Import Modules
import numpy as np
import pandas as pd

def sqtab(sqanti_out, ensg_to_gene, enst_to_trans, tabpath = "./sqanti_isoform_tab.tsv"):
    # Import Jurkat sqanti Output Text File
    cols = ['isoform', 'length', 'structural_category','associated_gene','associated_transcript','subcategory', 'FL', 'coding','ORF_length', 'CDS_start','CDS_end', 'predicted_NMD']
    data = pd.read_csv(sqanti_out, delimiter="\t", usecols = cols)
    data.columns = ['pb_acc', 'len', 'cat', 'gene','transcript', 'cat2', 'fl_cts','coding', 'orf_len', 'cds_st','cds_end', 'nmd']

    # Convert Structural Categories to Acronyms
    data.replace({"novel_not_in_catalog":"NNC","novel_in_catalog":"NIC","incomplete-splice_match":"ISM","full-splice_match":"FSM"}, inplace = True)

    # Filter out any cat that is not FSM, ISM, NNC or NIC
    f= ['FSM', 'ISM', 'NNC', "NIC"]
    fdata = data[data.cat.isin(f)]

    # Normalize fl_cts to cpm 
    sum = fdata['fl_cts'].sum(skipna=True)
    fdata['cpm'] = 1000000*fdata['fl_cts']/sum

    # Define 5'UTR and 3'UTR
    UTR5_len = fdata['cds_st'] - 1
    UTR3_len = fdata['len'] - fdata['cds_end'] + 1

    # Add 5'UTR and 3'UTR Information to Table
    fdata.insert(loc=len(fdata.columns), column = "5utr_len", value = UTR5_len)
    fdata.insert(loc=len(fdata.columns), column = "3utr_len", value = UTR3_len)

    ## Finding and Replacing Gene Information ##
    # Import Human Readable Gene Info and rename columns
    gen_name = pd.read_csv(ensg_to_gene, delimiter=r"\s+", header=None)
    gen_name.columns = ['A', 'B']

    # Make a Dictionary with Columns A and B
    gdict = pd.Series(gen_name.B.values,index=gen_name.A).to_dict()

    # Use Dictionary to Find Gene Names and Replace them with Human Readable Genes
    df = fdata[['gene']]
    fdata['gene'] = fdata['gene'].map(gdict).fillna(fdata['gene'])
    fdata.drop(fdata[fdata['gene'] == df['gene']].index, inplace=True)
  
    ## Finding and Replacing Transcript Information ##
    # Import Human Readable Transcript Info and rename columns
    trans = pd.read_csv(enst_to_trans, delimiter=r"\s+", header=None)
    trans.columns = ['A', 'B']

    # Make a Dictionary and Replace Transcript Names
    tdict = pd.Series(trans.B.values, index=trans.A).to_dict()
    fdata['transcript'] = fdata['transcript'].map(tdict).fillna(fdata['transcript'])

    # Write Dataframe as TSV File
    fdata.to_csv(tabpath, sep="\t", index= False, na_rep='0')

    print("Isoform Table from sqanti output has been prepared")

    # Sanity Check 
    #pdata.isin(['ENSG00000287807.1']).any()

    return fdata




