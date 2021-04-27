
MIN_ORF=30
ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results
python refine_orf_database.py \
--name jurkat_${MIN_ORF} \
--orfs $ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/filter_sqanti/src/jurkat_corrected.5degfilter.fasta
