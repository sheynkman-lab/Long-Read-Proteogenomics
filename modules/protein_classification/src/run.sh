PIPELINE_DIR=/Users/bj8th/Documents/Sheynkman-Lab/Data/21-05-05_Jurkat/jurkat

python protein_classification_add_meta.py \
--protein_classification  $PIPELINE_DIR/5p_utr/jurkat.sqanti_protein_classification_w_5utr_info.tsv \
--best_orf $PIPELINE_DIR/orf_calling/jurkat_best_orf.tsv \
--refined_meta $PIPELINE_DIR/refined_database/jurkat_orf_refined.tsv \
--name jurkat \
--dest_dir results/


python protein_classification.py \
--sqanti_protein results/jurkat.protein_classification_w_meta.tsv \
--name jurkat \
--dest_dir results/

