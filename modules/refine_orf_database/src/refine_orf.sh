MIN_ORF=20
ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results_chr1
python refine_orf_database.py \
--name jurkat_${MIN_ORF} \
--orfs $ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta

MIN_ORF=30
ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results_chr1
python refine_orf_database.py \
--name jurkat_${MIN_ORF} \
--orfs $ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta

MIN_ORF=40
ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results_chr1
python refine_orf_database.py \
--name jurkat_${MIN_ORF} \
--orfs $ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta

MIN_ORF=50
ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results_chr1
python refine_orf_database.py \
--name jurkat_${MIN_ORF} \
--orfs $ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta

MIN_ORF=30
ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results_chr1
python refine_orf_database.py \
--name jurkat_${MIN_ORF} \
--orfs $ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/filter_sqanti/src/jurkat_chr1.collapsed.fasta
