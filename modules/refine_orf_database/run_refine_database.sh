RESULTS=/Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/Long-Read-Proteogenomics/data/results/jurkat_gloria/results/results

python ./src/refine_orf_database.py \
  --name jurkat_1_cpm_filtered \
  --orfs $RESULTS/orf_calling/jurkat_best_orf_1_cpm.tsv \
  --pb_fasta $RESULTS/sqanti3-filtered/filtered_jurkat_corrected.fasta \
  --coding_score_cutoff 0.0 \
  