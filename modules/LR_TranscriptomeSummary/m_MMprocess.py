# module with functions related to processing of MetaMorpheus data

#### READ ME ####
## The function MMproc processes pacbio, genecode and uniport proteomics data from the MM output AllProteinGroups.tsv. 
## 
##
## MMproc(mm_out, db) = proc:
##                        inputs:
##                               mm_out = path of MM output AllProteins file to be processed
##                               db = 'pacbio', 'genecode' or 'uniprot'
##                               For db = genecode:
##                                      enst_to_gene file
##                               For db = uniprot:
##                                      uniprot_gene_to_genecode_gene file
##                        optional inputs:
##                               qvalue = For controlling FDR analyzed
##                                        DEFAULT: 0.01, 1% FDR
##                        outputs:
##                                proc = a dataframe with 'gene', 'psm', 'pep', 'upep', 'prot' and 'protsz'
##
## NOTE: Assumes that Genecode Gene mapping in MM output is not complete
## NOTE: Used Pacbio Gene Column 

## Import Modules ##
import numpy as np
import pandas as pd
from pathlib import Path
import csv

def MMproc(mm_out, db, trans_to_gene_file='./trans_to_gene.tsv', pbacc_to_gene_file = './uniprot_acc_to_gencode_gene.tsv',  qvalue = 0.01):
    # Check to see if mm_out exists
    pfile = Path(mm_out)
    if pfile.is_file()==False:
            print('The MM output does not exist. Please check your file path.')
    
    # Import MM output #
    if db == 'genecode':
        cols = ['Protein Accession', 'Number of Proteins in Group', 'Number of Peptides','Number of Unique Peptides',
        'Number of PSMs', 'Protein Decoy/Contaminant/Target', 'Protein QValue']
        col_names = ['pb_acc', 'prot', 'pep', 'upep', 'psm', 'dct', 'qval']
    else:
        cols = ['Protein Accession', 'Gene', 'Number of Proteins in Group', 'Number of Peptides','Number of Unique Peptides',
        'Number of PSMs', 'Protein Decoy/Contaminant/Target', 'Protein QValue']
        col_names = ['pb_acc','gene', 'prot', 'pep', 'upep', 'psm', 'dct', 'qval']
    
    data = pd.read_csv(mm_out, delimiter=r"\t", usecols = cols)
    data.columns = col_names 

    # Filter Data #
    fdata = data[data['qval'] <= qvalue]
    tdata = fdata[fdata['dct'] == 'T']

    # Map enst -> genes
    if db == 'genecode':
        # Make enst -> gene dictionary 
        trans_to_gene = pd.read_csv(trans_to_gene_file, sep='\t')
        trans_to_gene.columns = ['A', 'B']
        gdict = pd.Series(trans_to_gene.A.values, index = trans_to_gene.B).to_dict()


        # Extract pb_acc col and split each row 
        cut = tdata[['pb_acc']]
        split = cut.pb_acc.str.split('\||\.', expand=True)

        # Replace pb_acc -> gene_name 
        # TODO: There are some transcripts that don't map to genes. Ignoring those for now
        gen = split.apply(lambda x: x.map(gdict, na_action='ignore'))

    # Map uniprot_genes -> genecode_genes
    if db == 'uniprot':
        pbacc_to_gene = pd.read_csv(pbacc_to_gene_file, '\t')
        gdict = pd.Series(pbacc_to_gene.uniprot_gene.values, index =pbacc_to_gene.gencode_genes).to_dict()

        cut = tdata[['gene']]
        split = cut.gene.str.split('\||\.', expand=True)

        gen = split.apply(lambda x: x.map(gdict, na_action='ignore'))
    
    # Extract and split each gene from pacbio
    if db == 'pacbio':
        cut = tdata[['gene']]
        gen = cut.gene.str.split('\||\.', expand=True)
    
    # Remove duplicates
    tdata['gene'] = gen.stack().groupby(level=0).apply(lambda x: x.unique().tolist())

    # Extract genes and separate each unique gene name 
    p = tdata.gene.apply(pd.Series)

    # Insert psms, pep, upep and protein groups
    p.insert(0,'psm', tdata.psm.values)
    p.insert(0, 'pep', tdata.pep.values)
    p.insert(0, 'upep', tdata.upep.values)
    p.insert(0, 'prot', tdata.prot.values)

    # Sort by gene
    sort = p.melt(id_vars=['psm', 'pep', 'upep', 'prot'],value_name="gene").dropna().reset_index(drop=True).drop('variable',1)

    # Group by gene and sum psms, pep, upep and protein groups per gene
    if db == 'genecode':
        pro = sort.groupby(['gene']).agg(gc_psm=('psm','sum'), gc_pep=('pep', 'sum'), gc_upep=('upep', 'sum'), gc_prot=('prot', 'sum'), gc_protsz=('prot', 'mean')).reset_index().round(0)
    if db == 'uniprot':
        pro = sort.groupby(['gene']).agg(un_psm=('psm','sum'), un_pep=('pep', 'sum'), un_upep=('upep', 'sum'), un_prot=('prot', 'sum'), un_protsz=('prot', 'mean')).reset_index().round(0)
    if db == 'pacbio':
        pro = sort.groupby(['gene']).agg(pb_psm=('psm','sum'), pb_pep=('pep', 'sum'), pb_upep=('upep', 'sum'), pb_prot=('prot', 'sum'), pb_protsz=('prot', 'mean')).reset_index().round(0)

    pro['gene'] = pro['gene'].str.replace('-', '_')
    return pro 

