# run three times to compare databases pairwise

# gencode versus pacbio
python ./src/protein_groups_compare.py --pg_fileOne ./ProteinInfResults/GENCODE/AllQuantifiedProteinGroups.tsv --pg_fileTwo ./ProteinInfResults/PacBio/AllQuantifiedProteinGroups.tsv --mapping ./accession_map_gencode_uniprot_pacbio.tsv --output ./

# gencode versus uniprot
python ./src/protein_groups_compare.py --pg_fileOne <gencode protein groups file> --pg_fileTwo <uniprot protein groups> --mapping ./accession_map_gencode_uniprot_pacbio.tsv --output ./

# uniprot versus pacbio
python ./src/protein_groups_compare.py --pg_fileOne <uniprot protein groups file> --pg_fileTwo <pacbio protein groups> --mapping ./accession_map_gencode_uniprot_pacbio.tsv --output ./

