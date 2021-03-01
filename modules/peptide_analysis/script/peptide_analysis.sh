python ../src/peptide_analysis.py \
-gmap /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/reference_tables/gene_isoname.tsv \
-gc /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/metamorpheus/AllPeptides_Gencode.psmtsv \
-pb /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/refined_database/jurkat_orf_aggregated.fasta \
-sft /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/pacbio_6frm_gene_grouped/jurkat.6frame.fasta \
--pb_gene /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/transcriptome_summary/pb_gene.tsv \
--cpat_all_orfs /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/cpat/jurkat.ORF_prob.tsv \
--cpat_best_orf /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/cpat/jurkat.ORF_prob.best.tsv \
--cpat_orf_protein_fasta /Users/bj8th/Documents/Lab-for-Proteoform-Systems-Biology/LRPG-Manuscript/data/results/cpat/jurkat.ORF_seqs.fa \
-odir results/
