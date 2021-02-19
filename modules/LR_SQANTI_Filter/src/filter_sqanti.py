import pandas as pd 
import argparse
import logging
logging.basicConfig(filename='sqanti_filter.log', encoding='utf-8', level=logging.DEBUG) 

    """
    Filter SQANTI results based on several criteria 
    - protein coding only
    - percent A downstream
           perc_A_downstreamTTS : percent of genomic "A"s in the downstream 20 bp window. 
           If this number if high (say > 0.8), the 3' end site of this isoform is probably not reliable.
    - RTS stage
         RTS_stage: TRUE if one of the junctions could be a RT switching artifact.
    """

structural_categories = {
    'strict':['novel_not_in_catalog', 'novel_in_catalog',
        'incomplete-splice_match', 'full-splice_match', ]
    'all': : None
}


def string_to_boolean(string):
    """
    Converts string to boolean

    Parameters
    ----------
    string :str
    input string

    Returns
    ----------
    result : bool
    output boolean
    """
    if isinstance(string, bool):
        return str
    if string.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif string.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def filter_protein_coding(classification, protein_coding_filename):
    """
    Filter classification to only contain genes that are known to be protein-coding
    per Gencode
    
    Parameters
    ----------
    orfs : pandass DataFrame
        called ORFs
    protein_coding_filename : filename
        file of protien-coding genes. text file seperated by lines
    """
    logging.info("Filtering for only protein coding genes")
    with open(protein_coding_filename, 'r') as file:
        protein_coding_genes = file.read().splitlines()
    classification = classification[classification['gene'].isin(protein_coding_genes)]
    return classification


def filter_intra_polyA(classification):
    pass

def filter_rts_stage(classification):
    classification = classification[classification['RTS_stage'] == False]
    return classification


def get_filtered_classification(
    classification, 
    protein_coding_genes, 
    is_protein_coding_filtered, 
    is_intra_polyA_filtered,
    is_template_switching_filtered, 
    structural_categories_level ):
    
    classification = classification[~classification['associated_gene'].isna()]
    classification = classification[classification['associated_gene'].str.startswith("ENSG")]
    structural_categories_to_keep =['novel_not_in_catalog', 'novel_in_catalog',
        'incomplete-splice_match', 'full-splice_match', ]
    classification = classification[classification['structural_category'].isin(structural_categories_to_keep)]
    
    
    classification = classification[classification['perc_A_downstream_TTS'] <= percent_ployA_downstream]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", action="store", dest="name", help="sample name")
    parser.add_argument("--sqanti_classification", action='store', dest= 'classification_file',help='SQANTI Classification File')
    parser.add_argument("--sqanti_corrected_gtf", action="store", dest="corrected_gtf")
    parser.add_argument("--sqanti_corrected_fasta", action="store", dest="corrected_fasta")
    parser.add_argument("--filter_protein_coding", action="store", dest="filter_protein_coding", help"yes/no whether to keep only protein coding genes", default="yes")
    parser.add_argument("--filter_intra_polyA", action="store", dest="filter_intra_polyA", default="yes")
    parser.add_argument("--filter_template_switching", action="store", dest="filter_template_switching", default="yes")
    parser.add_argument("--protein_coding_genes", action="store", dest="protein_coding_genes", required=False)
    parser.add_argument("--percent_A_downstream_threshold", action="store", dest="percent_A_downstream_threshold", default=0.9,type=float)
    parser.add_argument("--structural_categories_level", action="store", dest="structural_categories_level", default="strict")
    results = parser.parse_args()

    is_protein_coding_filtered = string_to_boolean(results.filter_protein_coding)
    is_intra_polyA_filtered = string_to_boolean(results.filter_intra_polyA)
    is_template_switching_filtered = string_to_boolean(results.filter_template_switching)

    classification = pd.read_table(sqanti_classification_file)

    classification = get_filtered_classification()
    
    classification = classification[~classification['associated_gene'].isna()]
    classification = classification[classification['associated_gene'].str.startswith("ENSG")]

    if is_protein_coding_filtered:
        classification = filter_protein_coding(classification, results.protein_coding_genes)
    if is_intra_polyA_filtered:
        classification = filter_intra_polyA(classification)
    if is_template_switching_filtered:
        classification = filter_rts_stage(classification)

    





    






if __name__=="__main__":
    pass


# #%%
# canonical_chromisomes = [f"chr{i+1}" for i in range(22)]
# canonical_chromisomes = canonical_chromisomes + ['chrX','chrY','chrM']
# #%%
# import gtfparse
# import pandas as pd 

# gencode = gtfparse.parse_gtf("/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/Long-Read-Proteogenomics/data/input/gencode.v35.annotation.gtf")
# #%%

# with open("/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/Long-Read-Proteogenomics/data/input/gencode.v35.annotation.gtf") as gencode_in, open("/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/Long-Read-Proteogenomics/data/input/gencode.v35.annotation.canonical.gtf", "w") as gencode_out:
#     for line in gencode_in.readlines():
#         if line.startswith("#"):
#             gencode_out.write(line)
#         else:
#             chromisome = line.split("\t")[0].strip()
#             if chromisome in canonical_chromisomes:
#                 gencode_out.write(line)


# %%
