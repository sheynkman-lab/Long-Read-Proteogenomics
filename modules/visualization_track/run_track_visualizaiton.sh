
export RESULTS_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results
export REF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/input

nextflow run track_visualization.nf \
--name jurkat \
--include_transcript yes \
--sample_gtf $RESULTS_DIR/sqanti3-filtered/filtered_jurkat_corrected.gtf \
--refined_info $RESULTS_DIR/refined_database/jurkat_orf_refined.tsv \
--best_orf $RESULTS_DIR/orf_calling/jurkat_best_orf.tsv \
--pb_gene $RESULTS_DIR/transcriptome_summary/pb_gene.tsv \
--peptides $RESULTS_DIR/metamorpheus/pacbio/search_results/Task1SearchTask/AllPeptides.jurkat.psmtsv \
--refined_fasta $RESULTS_DIR/refined_database/jurkat_orf_refined.fasta \
--gencode_gtf $REF_DIR/gencode.v35.annotation.gtf 