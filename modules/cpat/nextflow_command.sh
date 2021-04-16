DATA_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/data
MIN_ORF=20
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta \
--min_orf $MIN_ORF \
--outdir results


MIN_ORF=30
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta \
--min_orf $MIN_ORF \
--outdir results

MIN_ORF=40
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta \
--min_orf $MIN_ORF \
--outdir results

MIN_ORF=50
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta \
--min_orf $MIN_ORF \
--outdir results

MIN_ORF=75
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/LRPG-Manuscript/data/results/sqanti3-filtered/filtered_jurkat_corrected.chr1.fasta \
--min_orf $MIN_ORF \
--outdir results


DATA_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/data
MIN_ORF=30
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/filter_sqanti/src/jurkat_chr1.collapsed.fasta \
--min_orf $MIN_ORF \
--outdir results_chr1