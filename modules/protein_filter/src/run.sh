# python filter_out_intergenic_and_cand_trunc.py \
# --protein_classification /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/protein_classification/input/jurkat_30_final.protein_classification.tsv \
# --gencode_gtf /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/input/gencode.v35.annotation.gtf \
# --protein_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/refined_database/jurkat_orf_refined.fasta \
# --name jurkat_30
NAME=jurkat
PIPELINE_DIR=/Users/bj8th/Documents/Sheynkman-Lab/Data/21-05-05_Jurkat/jurkat

python protein_filter.py \
--protein_classification /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/protein_classification/src/results/jurkat.protein_classification.tsv \
--gencode_gtf /Users/bj8th/Documents/Sheynkman-Lab/Data/Reference/gencode.v35.annotation.gtf \
--protein_fasta $PIPELINE_DIR/refined_database/${NAME}_orf_refined.fasta \
--sample_cds_gtf $PIPELINE_DIR/pacbio_cds/jurkat_with_cds.gtf \
--min_junctions_after_stop_codon 2 \
--name results/jurkat