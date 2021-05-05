echo "1_get_gc_exon_and_5utr_info"
python 1_get_gc_exon_and_5utr_info.py \
--gencode_gtf /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/input/gencode.v35.annotation.gtf \
--odir ./results

echo "2_classify_5utr_status"
python 2_classify_5utr_status.py \
--gencode_exons_bed results/gencode_exons_for_cds_containing_ensts.bed \
--gencode_exons_chain results/gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
--sample_cds_gtf /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/pacbio_cds/jurkat_with_cds.gtf \
--odir ./results 

echo "3_merge_5utr_info_to_pclass_table"
python 3_merge_5utr_info_to_pclass_table.py \
--name jurkat \
--utr_info results/pb_5utr_categories.tsv \
--sqanti_protein_classification /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/protein_classification/input/jurkat_30.sqanti_protein_classification.tsv \
--odir ./results