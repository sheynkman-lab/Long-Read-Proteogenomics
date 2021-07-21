#!/bin/bash
#-----------------------------------------------------------------------------------------------------------
#
# Convenience documentation and script for downloading and when and if necessary reconstruct the files for
# proper execution of the nextflow pipeline.
#
# Data will be downloaded into the Long-Read-Proteogenomics/data directory.
# As part of the Sheynkman lab's long read proteogenomics pipeline.
#
# It is assumed that the https://github.com/sheynkman-lab/Long-Read-Proteogenomics.git has been cloned
#
# Please run from the Long-Read-Proteogenomics/data directory
#
# These data previously were stored in Amazon S3 buckets.
# And now reside in Zenodo
#
# The main reasons why ZENODO vs AWS S3 is a good idea:
#    1. Data versioning (number 1 important reason),
#       In S3 data can be overwritten for the same path at any point possibly breaking the pipeline.
#    2. Cost, removing data from S3. These data are tiny but the principle stays: The less storage the better
#    3. Access, Academic users familiarity/accessibility - Most reviewers, readers of the pipeline and paper
#       will know ZENODO and will be able to use the data, AWS has an entry barrier for this group of people
#
#-----------------------------------------------------------------------------------------------------------
echo "downloading from Zenodo record 5076056\n\n"

wget https://zenodo.org/record/5076056/files/gencode.v35.annotation.canonical.gtf
wget https://zenodo.org/record/5076056/files/gencode.v35.pc_transcripts.fa.gz
wget https://zenodo.org/record/5076056/files/gencode.v35.pc_translations.fa.gz
wget https://zenodo.org/record/5076056/files/GRCh38.primary_assembly.genome.canonical.fa.gz
wget https://zenodo.org/record/5076056/files/Human_Hexamer.tsv
wget https://zenodo.org/record/5076056/files/Human_logitModel.RData
wget https://zenodo.org/record/5076056/files/jurkat_classification.txt
wget https://zenodo.org/record/5076056/files/jurkat_corrected.fasta.gz
wget https://zenodo.org/record/5076056/files/jurkat_corrected.gtf
wget https://zenodo.org/record/5076056/files/jurkat_gene_kallisto.tsv
wget https://zenodo.org/record/5076056/files/jurkat_merged.ccs.bam
wget https://zenodo.org/record/5076056/files/jurkat_r1.fastq.gz
wget https://zenodo.org/record/5076056/files/jurkat_r2.fastq.gz
wget https://zenodo.org/record/5076056/files/kallist_table_rdeplete_jurkat.tsv
wget https://zenodo.org/record/5076056/files/mass_spec.tar.gz
wget https://zenodo.org/record/5076056/files/NEB_primers.fasta
wget https://zenodo.org/record/5076056/files/star_genome.tar.gz
wget https://zenodo.org/record/5076056/files/Task1SearchTaskconfig_orf.toml
wget https://zenodo.org/record/5076056/files/uniprot_reviewed_canonical_and_isoform.fasta.gz

echo "done downloading ... now reconstructing folders untarring and  unzipping mass_spec.tar.gz and star_genome.tar.gz\n\n"

tar xvzf mass_spec.tar.gz
tar xvzf star_genome.tar.gz

echo "done reconstructing"
