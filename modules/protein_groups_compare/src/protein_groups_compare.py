#!/usr/bin/env python3

# determine the relationship between protein groups identified using gencode, uniprot, and/or pacbio databases

# %%

import pandas as pd
import numpy as np
from collections import defaultdict
import argparse


# function to remove decoy protein groups, protein groups above a 1% FDR cutoff, and extraneous columns 
def tsv_filter(tsv, col_to_keep=['Protein Accession', 'Gene', 'Unique Peptides','Shared Peptides','Sequence Coverage Fraction','Number of PSMs','Protein QValue']):
    """
    Filters the dataframe for FDR, target proteins and removes columns which are not pertainant
    
    Parameters
    -----------
    tsv : pandas DataFrame
        loaded AllProteinGroups.tsv file from MetaMorpheus
    col_to_keep: list of str
        columns which will remain in the dataframe

    Returns
    -----------
    tsvfilter: pandas DataFrame    
    """   
    tsvfilter = tsv
    tsvfilter = tsvfilter.loc[(tsvfilter['Protein QValue'] <= 0.01)] #filter for 1%FDR
    tsvfilter = tsvfilter.loc[(tsvfilter['Protein Decoy/Contaminant/Target'] == "T")]
    for col in tsv.columns:
        if col not in col_to_keep:
            del tsvfilter[col]

    return tsvfilter

def read_filter_proteinGroups(filePath):
    """
    Reads the AllProteinGroups.tsv file in, filters the data and columns.
    Keep only the Target protein groups below 1% FDR, and only columns that contain valuable information to carry forward

    Parameters
    ------------
    filePath: str
        location of the AllProteinGroups.tsv file

    Returns
    -----------
    proteinGroups_filtered : pandas DataFrame
    """
    proteinGroups = pd.read_table(filePath, index_col=False)
    proteinGroups_filtered = tsv_filter(proteinGroups)
    return proteinGroups_filtered


def format_accession_mapping_table(conversionKey):
    """
    Reformats the accession mapping file to be manipulated going forward

    Parameters
    -----------
    conversionKey: pandas DataFrame
        loaded accession mapping file

    Return
    -----------
    imap: pandas DataFrame 
    """
    isomap = conversionKey
    isomap['un'] = isomap['uniprot_acc']
    isomap['pb'] = isomap['pacbio_acc']
    isomap['gc'] = isomap['gencode_acc']
    imap = isomap[['pb', 'un', 'gc']].dropna(how="all")
    imap['idx'] = np.arange(len(imap))

    return imap


def return_list_with_nans_filtered_out(input_list):
    '''
    return_list_with_nans_filtered_out - function to remove NAs from a list
    #Input:
    input_list - list - list to remove NAs from
    #Returns:
    filtered_list - list - list without NAs
    '''
    filtered_list = []
    for item in input_list:
        if pd.isnull(item): continue
        filtered_list.append(item)
    return filtered_list

def accession_mapping_setup (filePath):
    """
    Reads in the accession mapping key from .tsv and creates conversion dictionaries for a PacBio and UniProt starting point

    Parameters
    ------------
    filePath: str
        location of the accession mapping file

    Returns
    ------------
    conversionDics = list of dictionaries
    """
    accession_conversion_key = pd.read_table(filePath)
    #format conversion key
    
    imap = format_accession_mapping_table (accession_conversion_key)

    # create pb to ref accession map
    pb_ref = defaultdict() # pb_acc -> [pb, un, gc] (no entry if no accession)
    un_ref = defaultdict() # un_acc -> [pb, un, gc] (no entry if no accession)

    for i, row in imap.iterrows():
        pb, un, gc, idx = row
        acc_list = return_list_with_nans_filtered_out([pb, un, gc])
        pb_ref[pb] = acc_list[-1]
        un_ref[un] = acc_list[-1]

    dictionaryList = []
    dictionaryList.append(pb_ref)
    dictionaryList.append(un_ref)

    return dictionaryList

#identify databases as pacbio (PacBio), uniprot (UniProt), or gencode (GENCODE)
def grp_is_derived_from_pacbio(grp):
    '''
    grp_is_derived_from_pacbio - function to determine whether an accession number comes from pac-bio
    #Input:
    grp - str - the accession number
    #Returns:
    boolean - True if from pac-bio otherwise False
    '''
    peek_acc = grp 
    if peek_acc.startswith('PB.'):
        return True
    else:
        return False

def grp_is_derived_from_gencode(grp):
    '''
    grp_is_derived_from_gencode - function to determine whether an accession number comes from gencode
    #Input:
    grp - str - the accession number
    #Returns:
    boolean - True if from gencode otherwise False
    '''
    peek = grp 
    if '-' not in peek:
        if peek.startswith('ENSP'):
            return True
        else:            
            return False 
    else:
        iso_num = int(peek.split('-')[-1].split('|')[0])
        if iso_num > 100:
            return True
        else:
            return False

def determine_source_database_for_grps(unassigned_grps):
    '''
    determine_source_database_for_grps - function to list what database an accession number comes from (gencode: GENCODE, uniprot: UniProt, pac-bio: PacBio)
    #Input:
    unassigned_grps - list - the accession numbers
    #Returns:
    grps - list - list identifying where each accession number came from
    '''
    grps = [] 
    for grp in unassigned_grps:
        if grp_is_derived_from_pacbio(grp):
            grps.append('PacBio') #grps.append(['PacBio':grp])
        elif grp_is_derived_from_gencode(grp):
            grps.append('GENCODE') #grps.append(['GENCODE':grp])
        else:
            grps.append('UniProt') #grps.append(['UniProt':grp])
    
    return grps

#select most common string as assignment
def most_frequent(search_list):
    '''
    most_frequent - function to determine the most frequently used string in a list
    #Input:
    search_list - list - list of strings
    #Returns:
    most frequently appearing string
    '''
    return max(set(search_list), key = search_list.count)

#add accession number conversion column 
def translate_Accessions (pg, db_id, conversionDictionaries):
    pb_ref = conversionDictionaries[0]
    un_ref = conversionDictionaries[1]
    all_trans = []
    for proteinGroup in pg['Protein Accession']:
        groupTrans = set()
        if '|' not in proteinGroup:
            if db_id == 'UniProt':
                if proteinGroup in un_ref:
                    groupTrans.add(un_ref[proteinGroup])
                else:
                    groupTrans.add(proteinGroup)
            if db_id == 'PacBio':
                if proteinGroup in pb_ref:
                    groupTrans.add(pb_ref[proteinGroup])
                else:
                    groupTrans.add(proteinGroup)
            if db_id == 'GENCODE':
                groupTrans.add(proteinGroup)
        else:
            split_group = [str(proteinGroup).split('|')]
            for accession in split_group:
                for split_accession in accession:
                    if db_id == 'UniProt':
                        if split_accession in un_ref:
                            groupTrans.add(un_ref[split_accession])
                        else:
                            groupTrans.add(split_accession)
                    if db_id == 'PacBio':
                        if split_accession in pb_ref:
                            groupTrans.add(pb_ref[split_accession])
                        else:
                            groupTrans.add(split_accession)
                    if db_id == 'GENCODE':
                        groupTrans.add(split_accession)
        all_trans.append(groupTrans)
    pg['Protein Accession Translation'] = all_trans



#compare the protein groups for the two different analyses
def accession_compare(df1,df2,db1_id, db2_id, output):
    #columns for the cumulative dataframes that will be exported (basic info for each protein group with the database of origin appended)
    cols = ['Protein Accession {}'.format(db1_id), 'Gene {}'.format(db1_id), 'Unique Peptides {}'.format(db1_id),'Shared Peptides {}'.format(db1_id),'Sequence Coverage Fraction {}'.format(db1_id),'Number of PSMs {}'.format(db1_id), 'Protein QValue {}'.format(db1_id), 'Protein Accession Translation {}'.format(db1_id), 'Protein Accession {}'.format(db2_id), 'Gene {}'.format(db2_id), 'Unique Peptides {}'.format(db2_id),'Shared Peptides {}'.format(db2_id),'Sequence Coverage Fraction {}'.format(db2_id),'Number of PSMs {}'.format(db2_id), 'Protein QValue {}'.format(db2_id), 'Protein Accession Translation {}'.format(db2_id)]
    #columns for the dataframes that are for protein groups unique to a single database
    cols_1=['Protein Accession {}'.format(db1_id), 'Gene {}'.format(db1_id), 'Unique Peptides {}'.format(db1_id),'Shared Peptides {}'.format(db1_id),'Sequence Coverage Fraction {}'.format(db1_id),'Number of PSMs {}'.format(db1_id), 'Protein QValue {}'.format(db1_id), 'Protein Accession Translation {}'.format(db1_id)]
    cols_2=['Protein Accession {}'.format(db2_id), 'Gene {}'.format(db2_id), 'Unique Peptides {}'.format(db2_id),'Shared Peptides {}'.format(db2_id),'Sequence Coverage Fraction {}'.format(db2_id),'Number of PSMs {}'.format(db2_id), 'Protein QValue {}'.format(db2_id), 'Protein Accession Translation {}'.format(db2_id)]
    # setting up empty dataframes that will be filled with categories of protein group matchings
    No_shared_accessions_db1=pd.DataFrame(columns=cols_1)
    No_shared_accessions_db2=pd.DataFrame(columns=cols_2)
    Exactly_one_match=pd.DataFrame(columns=cols)
    Simpler_in_one=pd.DataFrame(columns=cols)
    Simpler_in_two=pd.DataFrame(columns=cols)
    Partially_overlapping=pd.DataFrame(columns=cols)

    #re-label the dataframe columms with database of origin, and remove T/D/C column
    df1.columns=['Protein Accession {}'.format(db1_id), 'Gene {}'.format(db1_id), 'Unique Peptides {}'.format(db1_id),'Shared Peptides {}'.format(db1_id),'Sequence Coverage Fraction {}'.format(db1_id),'Number of PSMs {}'.format(db1_id), 'Protein QValue {}'.format(db1_id), 'Protein Accession Translation {}'.format(db1_id)]
    df2.columns=['Protein Accession {}'.format(db2_id), 'Gene {}'.format(db2_id), 'Unique Peptides {}'.format(db2_id),'Shared Peptides {}'.format(db2_id),'Sequence Coverage Fraction {}'.format(db2_id),'Number of PSMs {}'.format(db2_id), 'Protein QValue {}'.format(db2_id), 'Protein Accession Translation {}'.format(db2_id)]
    
    df1 = df1.dropna(how="all")
    df2 = df2.dropna(how="all")
    
    #list of sets (accessions in a protein group) that will be compared
    list_set_db1 = []
    list_set_db2 = []    
    
    #fill the lists of sets
    for accession in df1['Protein Accession Translation {}'.format(db1_id)]:     
        list_set_db1.append(accession)
    
    for accession in df2['Protein Accession Translation {}'.format(db2_id)]:
        list_set_db2.append(accession) 

    #list of sets from db2 that have found a match (helps us later un figure out which groups are unique to db2)      
    list_set_db2_matches =[]
    for group1 in list_set_db1:        
        
        distinct = True      
        for group2 in list_set_db2:
            if len(group1.intersection(group2))>0:
                distinct= False
                if len(group1.intersection(group2)) == len(group1.union(group2)):
                    # add to exact match
                    data1 = df1.loc[df1['Protein Accession Translation {}'.format(db1_id)] == group1]                    
                    data2 = df2.loc[df2['Protein Accession Translation {}'.format(db2_id)] == group2]
                    Exactly_one_match=Exactly_one_match.append(pd.DataFrame([[data1.iat[0, 0], data1.iat[0,1], data1.iat[0, 2], data1.iat[0, 3], data1.iat[0, 4], data1.iat[0, 5], data1.iat[0, 6], data1.iat[0, 7],data2.iat[0, 0], data2.iat[0,1], data2.iat[0, 2], data2.iat[0, 3], data2.iat[0, 4], data2.iat[0, 5], data2.iat[0, 6], data2.iat[0, 7]]], columns=cols))
                    list_set_db2_matches.append(group2)                    
                    break
                elif group1.issubset(group2):
                    # is simpler in group1
                    data1 = df1.loc[df1['Protein Accession Translation {}'.format(db1_id)] == group1]
                    data2 = df2.loc[df2['Protein Accession Translation {}'.format(db2_id)] == group2]
                    Simpler_in_one=Simpler_in_one.append(pd.DataFrame([[data1.iat[0, 0], data1.iat[0,1], data1.iat[0, 2], data1.iat[0, 3], data1.iat[0, 4], data1.iat[0, 5], data1.iat[0, 6], data1.iat[0, 7],data2.iat[0, 0], data2.iat[0,1], data2.iat[0, 2], data2.iat[0, 3], data2.iat[0, 4], data2.iat[0, 5], data2.iat[0, 6], data2.iat[0, 7]]], columns=cols))
                    list_set_db2_matches.append(group2)                    
                    break
                elif group1.issuperset(group2):
                    # is simpler in group2
                    data1 = df1.loc[df1['Protein Accession Translation {}'.format(db1_id)] == group1]
                    data2 = df2.loc[df2['Protein Accession Translation {}'.format(db2_id)] == group2]
                    Simpler_in_two=Simpler_in_two.append(pd.DataFrame([[data1.iat[0, 0], data1.iat[0,1], data1.iat[0, 2], data1.iat[0, 3], data1.iat[0, 4], data1.iat[0, 5], data1.iat[0, 6], data1.iat[0, 7],data2.iat[0, 0], data2.iat[0,1], data2.iat[0, 2], data2.iat[0, 3], data2.iat[0, 4], data2.iat[0, 5], data2.iat[0, 6], data2.iat[0, 7]]], columns=cols))
                    list_set_db2_matches.append(group2)                    
                else:
                    # overlap
                    data1 = df1.loc[df1['Protein Accession Translation {}'.format(db1_id)] == group1]
                    data2 = df2.loc[df2['Protein Accession Translation {}'.format(db2_id)] == group2]
                    Partially_overlapping=Partially_overlapping.append(pd.DataFrame([[data1.iat[0, 0], data1.iat[0,1], data1.iat[0, 2], data1.iat[0, 3], data1.iat[0, 4], data1.iat[0, 5], data1.iat[0, 6], data1.iat[0, 7],data2.iat[0, 0], data2.iat[0,1], data2.iat[0, 2], data2.iat[0, 3], data2.iat[0, 4], data2.iat[0, 5], data2.iat[0, 6], data2.iat[0, 7]]], columns=cols))
                    list_set_db2_matches.append(group2)                    
        if distinct == True:
            #completely distinct
            data1 = df1.loc[df1['Protein Accession Translation {}'.format(db1_id)] == group1]                   
            No_shared_accessions_db1=No_shared_accessions_db1.append([data1])            

    list_set_unique_db2 = np.setdiff1d(list_set_db2, list_set_db2_matches)
    for group in list_set_unique_db2:                    
        data = df2.loc[df2['Protein Accession Translation {}'.format(db2_id)] == group]
        No_shared_accessions_db2=No_shared_accessions_db2.append([data])

    #this is the one output file of the whole module (will not be used as input for any other module) 
    fileName = '/ProteinInference_{}_{}_comparisons.xlsx'.format(db1_id, db2_id)
    outputPath = output + fileName  
    writer = pd.ExcelWriter(outputPath, engine='xlsxwriter')
    Exactly_one_match.to_excel(writer, sheet_name='ProteinGroup_ExactMatches',index = False)
    Simpler_in_one.to_excel(writer, sheet_name='ProteinGroup_SimplerIn_{}'.format(db1_id),index = False)
    Simpler_in_two.to_excel(writer, sheet_name='ProteinGroup_SimplerIn_{}'.format(db2_id),index = False)
    Partially_overlapping.to_excel(writer, sheet_name='ProteinGroup_PartiallyOverlap',index = False)
    No_shared_accessions_db1.to_excel(writer, sheet_name='ProteinGroup_DistinctTo_{}'.format(db1_id),index = False)
    No_shared_accessions_db2.to_excel(writer, sheet_name='ProteinGroup_DistinctTo_{}'.format(db2_id),index = False)
    writer.save()      

def main():
    parser= argparse.ArgumentParser(description ='Compare Protein Inference Results from 2 different MetaMorpheus Searches')
    parser.add_argument('--pg_fileOne', '-pg_1', action = 'store', dest ='proteinGroups_db1', help = 'File location of the AllProteinGroups.tsv file from the first MetaMorpheus search')
    parser.add_argument('--pg_fileTwo', '-pg_2', action = 'store', dest = 'proteinGroups_db2', help = 'File location of the AllProteinGroups.tsv file from the second MetaMorpheus search')
    parser.add_argument('--mapping', '-m', action = 'store', dest = 'accession_conversion_key', help = 'File location of the Accession mapping file for the converting accessions for comparison purposes')
    parser.add_argument('--output', '-o', action = 'store', dest = 'output', help = 'Output file location')
    scriptInput = parser.parse_args()

    proteinGroups_db1_filtered = read_filter_proteinGroups(scriptInput.proteinGroups_db1)
    proteinGroups_db2_filtered = read_filter_proteinGroups(scriptInput.proteinGroups_db2)

    conversionDictionaries = accession_mapping_setup(scriptInput.accession_conversion_key)
   
    #Determine the database of origin for the  forst 5 protein groups, and pick the most common (prevents issues with potential overlap between UniProt and GENCODE criteria)
    db1_id = most_frequent(determine_source_database_for_grps(proteinGroups_db1_filtered.head(5)['Protein Accession']))
    db2_id = most_frequent(determine_source_database_for_grps(proteinGroups_db2_filtered.head(5)['Protein Accession']))

    proteinGroups_db1_filtered = proteinGroups_db1_filtered.dropna(how = "all")
    proteinGroups_db2_filtered = proteinGroups_db2_filtered.dropna(how ="all")

    #Adds column to protein group DataFrame containing a set of translated protein groups for each grouping
    translate_Accessions(proteinGroups_db1_filtered, db1_id, conversionDictionaries)
    translate_Accessions(proteinGroups_db2_filtered, db2_id, conversionDictionaries)

    accession_compare(proteinGroups_db1_filtered, proteinGroups_db2_filtered, db1_id, db2_id, scriptInput.output)
 
if __name__ == "__main__":
    main()

