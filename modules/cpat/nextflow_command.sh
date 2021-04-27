DATA_DIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/data

MIN_ORF=30
nextflow cpat.nf \
--name jurkat_${MIN_ORF} \
--hexamer $DATA_DIR/input/Human_Hexamer.tsv \
--logit_model $DATA_DIR/input/Human_logitModel.RData \
--sample_fasta /Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics/modules/filter_sqanti/src/jurkat_corrected.5degfilter.fasta \
--min_orf $MIN_ORF \
--outdir results