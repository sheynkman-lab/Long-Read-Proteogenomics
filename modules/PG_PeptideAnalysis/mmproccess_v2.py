"""
This module prepares a summary table from the output of mass spec searches

    Inputs:
    ----------------------------------------------------------------------
    1. mmout: AllProteins file from MM output 
    
    For database = gencode:
        2. enst -> gene file 
    
    For database = uniprot:
        2. Uniprot gene -> Gencode gene file 
    -----------------------------------------------------------------------

    Optional Inputs:
    -----------------------------------------------------------------------
    1. qvalue: Control % FDR, default = 1% FDR
    ------------------------------------------------------------------------

    Outputs:
    -----------------------------------------------------------------------
    1. proc: a table of a gene and its associated psm, pep, unique pep and protein groups
    -----------------------------------------------------------------------
"""

