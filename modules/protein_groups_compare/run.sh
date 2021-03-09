# run three times to compare databases pairwise
RESULTS=/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/Long-Read-Proteogenomics/data/results/jurkat_gloria/results/results
# gencode versus pacbio
echo 'PROTEIN GROUPS COMPARE: GENCODE PACBIO'
python ./src/protein_groups_compare.py --pg_fileOne $RESULTS/metamorpheus/gencode/search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv --pg_fileTwo $RESULTS/metamorpheus/pacbio/search_results/Task1SearchTask/AllQuantifiedProteinGroups.jurkat.tsv --mapping $RESULTS/accession_mapping/accession_map_gencode_uniprot_pacbio.tsv --output $RESULTS/protein_groups_compare

# # gencode versus uniprot
# echo 'PROTEIN GROUPS COMPARE: GENCODE UNIPROT'
python ./src/protein_groups_compare.py --pg_fileOne $RESULTS/metamorpheus/uniprot/search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv --pg_fileTwo $RESULTS/metamorpheus/gencode/search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv --mapping $RESULTS/accession_mapping/accession_map_gencode_uniprot_pacbio.tsv --output $RESULTS/protein_groups_compare

# # uniprot versus pacbio
# echo 'PROTEIN GROUPS COMPARE: UNIPROT PACBIO'
python ./src/protein_groups_compare.py --pg_fileOne $RESULTS/metamorpheus/uniprot/search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv --pg_fileTwo $RESULTS/metamorpheus/pacbio/search_results/Task1SearchTask/AllQuantifiedProteinGroups.jurkat.tsv --mapping $RESULTS/accession_mapping/accession_map_gencode_uniprot_pacbio.tsv --output $RESULTS/protein_groups_compare

