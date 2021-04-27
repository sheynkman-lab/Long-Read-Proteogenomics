RESULTS_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results
REF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/input
BEST_ORF_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/orf_calling/src/results
REFINED_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/refine_orf_database/src/results

# for MIN_ORF in 20 30 40 50 75
# do
MIN_ORF=30
nextflow run track_visualization.nf \
--name jurkat_${MIN_ORF} \
--include_transcript yes \
--sample_gtf /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/filter_sqanti/src/jurkat_corrected.5degfilter.gff \
--refined_info $REFINED_DIR/jurkat_${MIN_ORF}_orf_refined.tsv \
--best_orf $BEST_ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
--pb_gene $RESULTS_DIR/transcriptome_summary/pb_gene.tsv \
--refined_fasta $REFINED_DIR/jurkat_${MIN_ORF}_orf_refined.fasta \
--peptides $RESULTS_DIR/metamorpheus/pacbio/search_results/Task1SearchTask/AllPeptides.jurkat.psmtsv \
--gencode_gtf $REF_DIR/gencode.v35.annotation.gtf \
--outdir ./results
# done




# # for MIN_ORF in 20 30 40 50 75
# # do
# MIN_ORF=30
# nextflow run track_visualization.nf \
# --name jurkat_${MIN_ORF} \
# --include_transcript yes \
# --sample_gtf $RESULTS_DIR/sqanti3-filtered/filtered_jurkat_corrected.chr1.gtf \
# --refined_info $REFINED_DIR/jurkat_${MIN_ORF}_orf_refined.tsv \
# --best_orf $BEST_ORF_DIR/jurkat_${MIN_ORF}_best_orf.tsv \
# --pb_gene $RESULTS_DIR/transcriptome_summary/pb_gene.tsv \
# --refined_fasta $REFINED_DIR/jurkat_${MIN_ORF}_orf_refined.fasta \
# --peptides $RESULTS_DIR/metamorpheus/pacbio/search_results/Task1SearchTask/AllPeptides.jurkat.psmtsv \
# --gencode_gtf $REF_DIR/gencode.v35.annotation.chr1.gtf \
# --outdir ./results_chr1 
# # done


