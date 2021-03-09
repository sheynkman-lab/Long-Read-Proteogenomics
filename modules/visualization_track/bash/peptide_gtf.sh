RESULTS=../../../data/results/jurkat_gloria/results/results
python ../src/make_peptide_gtf_file.py \
--name jurkat \
--sample_gtf $RESULTS/pacbio_cds/jurkat_with_cds.gtf \
--peptides $RESULTS/metamorpheus/pacbio/search_results/Task1SearchTask/AllPeptides.jurkat.psmtsv \
--pb_gene $RESULTS/transcriptome_summary/pb_gene.tsv \
--refined_fasta $RESULTS/refined_database/jurkat_orf_refined.fasta