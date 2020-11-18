
import argparse

import pandas as pd


def process_metamorpheus_file(filename):
    """
    Process MetaMorpheus AllPeptides file 
    Keeps only target peptieds with FDR <= 0.01
    Explodes 'Protein Accession' column into individual rows

    Parameters
    ----------
    filename : str
        MetaMorpheus AllPeptides file location
    
    Returns
    --------
    original : pandas DataFrame
        exploded AllPeptieds dataframe
    """
    original = pd.read_csv(filename, sep = '\t')
    below_fdr = original['QValue'] <= 0.01
    target = original['Decoy/Contaminant/Target'] == 'T'
    original = original[below_fdr & target]
    original['Protein Accession'] = original['Protein Accession'].str.split('|')
    original = original.explode('Protein Accession')
    return original


def find_best_protein(remaining, original, protein_column):
    """
    Finds best protein in remaining graph based on number
    of peptides a protein matches to (# of edges a protein has). Choose max
    If there are proteins with an equal number of maximum matched peptides
    then those proteins compared to the original graph for number of edges 
    a protein has
    """
    sizes = remaining.groupby(protein_column).size().reset_index().rename(columns = {0:'size'})
    max_size = sizes['size'].max()
    best_proteins = list(sizes[sizes['size'] == max_size ][protein_column])
    if len(best_proteins) > 1:
        subset = original[original[protein_column].isin(best_proteins)]
        subsizes = subset.groupby(protein_column).size().reset_index().rename(columns = {0:'size'})
        max_subsize = subsizes['size'].max()
        best_proteins = list(subsizes[subsizes['size'] == max_subsize][protein_column])
    return best_proteins[0]



def greedy_inference(original, protein_column = 'Protein Accession', peptide_column = 'Base Sequence'):
    """
    Greedy protein inference algorithm for matching peptids to corresponding proteins

    Notaion:
    G : original graph
    Gi : inferred graph
    Gr : remaining graph
    Gd: dropped graph
    p : greedily selcted protein
    s : peptides connected to p


    Select peptides in G that only match to single protein
    Add proteins corresponding to peptides and all attached peptides to Gi
    Remove said proteins from  Gr
    While Gr has edges connected proteins and peptides
        Greedily select best protein p
        Add p and connected peptides Gi
        Add peptide-protein edges where protein is not p and peptide is in s in Gd
        Remove edgees where peptides is in s from Gr
    
    Remake Gi and make Gd
        Gi remade to contain all protein-peptide edges that connect to an inferred protein
        Gd made to contain all protein-peptide edges that do not connect to an inferred protein

    Parameters
    ---------
    original : pandas DataFrame
        original peptide-protien graph
    protein_column : str
        column associated with protein accession
    peptide_column : str
        column associated with peptide

    Returns
    --------
    inferred: pandas DataFrame
        Gi, subgraph of G of proteins and their associated peptides
    dropped: pandas DataFrame
        Gd, subgraph of G of proteins and their associated peptides

    """

    dropped = pd.DataFrame(columns = original.columns)
    peptide_sizes = original.groupby(peptide_column).size().reset_index().rename(columns = {0:'size'})
    single = peptide_sizes[peptide_sizes['size'] == 1]
    inferred = original[original[peptide_column].isin(single[peptide_column])]
    inferred = original[original[protein_column].isin(inferred[protein_column])]
    remaining = original[~original[protein_column].isin(inferred[protein_column])]
    inferred = [ inferred ]

    while len(remaining) > 0:
        best_protein = find_best_protein(remaining, original, protein_column)
        matches = remaining[remaining[protein_column] == best_protein]
        tmp_peptides = list(matches[peptide_column])
        inferred.append(matches)

        is_matched_peptide = remaining[peptide_column].isin(tmp_peptides)
        is_best_protein = remaining[protein_column] == best_protein
        
        remaining = remaining[~is_matched_peptide]
        
    inferred = pd.concat(inferred)

    inferred_proteins = inferred[protein_column].unique()
    inferred = original[original[protein_column].isin(inferred_proteins)]
    dropped = original[~original[protein_column].isin(inferred_proteins)]
    
    return inferred, dropped
 

def main():
    parser = argparse.ArgumentParser(description='Parse Greedy Protein Inference Arguments')
    parser.add_argument('-i','--ifile', action='store',dest='ifile', help = 'input peptide-protein file')
    parser.add_argument('-o','--odir', action = 'store', dest='odir', help = 'output directory for results')
    results = parser.parse_args()

    original = process_metamorpheus_file(results.ifile)
    inferred, dropped = greedy_inference(original)
    inferred.to_csv(f'{results.odir}/inferred_proteins.tsv', sep = '\t', index = False)
    dropped.to_csv(f'{results.odir}/dropped_proteins.tsv', sep = '\t', index=False)

if __name__ == "__main__":
    main()

        



