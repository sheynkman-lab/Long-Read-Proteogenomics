PIPELINE_DIR=/Users/bj8th/Documents/Sheynkman-Lab/Data/21-05-05_Jurkat/jurkat

python protein_database_aggregation.py \
--protein_classification /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/protein_filter/src/results/jurkat.classification_filtered.tsv \
--gene_lens $PIPELINE_DIR/reference_tables/gene_lens.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/protein_filter/src/results/jurkat.filtered_protein.fasta \
--gc_fasta $PIPELINE_DIR/gencode_db/gencode_protein.fasta \
--refined_info $PIPELINE_DIR/refined_database/jurkat_orf_refined.tsv \
--pb_cds_gtf $PIPELINE_DIR/pacbio_cds/jurkat_with_cds.gtf \
--name results/jurkat \
--lower_kb 1 \
--upper_kb 4 \
--lower_cpm 3 \