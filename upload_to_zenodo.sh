#!/bin/bash
#-----------------------------------------------------------------------------------------------------------
#
# Convenience documentation and script for uploading data files from
# Long-Read-Proteogenomics/data directory. As part of the Sheynkman lab's long read proteogenomics pipeline.
#
# These data previously were stored in Amazon S3 buckets.
# This upload puts the data in the Zenodo repository for ease of access
#
# The main reasons why ZENODO vs AWS S3 is a good idea:
#    1. Data versioning (number 1 important reason),
#       In S3 data can be overwritten for the same path at any point possibly breaking the pipeline.
#    2. Cost, removing data from S3. These data are tiny but the principle stays: The less storage the better
#    3. Access, Academic users familiarity/accessibility - Most reviewers, readers of the pipeline and paper
#       will know ZENODO and will be able to use the data, AWS has an entry barrier for this group of people
#
#-----------------------------------------------------------------------------------------------------------

cd ../zenodo-upload
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/gencode.v35.annotation.canonical.gtf
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/gencode.v35.pc_transcripts.fa.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/gencode.v35.pc_translations.fa.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/GRCh38.primary_assembly.genome.canonical.fa.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/Human_Hexamer.tsv
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/Human_logitModel.RData
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_classification.txt
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_corrected.fasta.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_corrected.gtf
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_gene_kallisto.tsv
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_merged.ccs.bam
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_r1.fastq.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/jurkat_r2.fastq.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/kallist_table_rdeplete_jurkat.tsv
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/mass_spec.tar.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/NEB_primers.fasta
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/star_genome.tar.gz
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/Task1SearchTaskconfig_orf.toml
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/uniprot_reviewed_canonical_and_isoform.fasta.gz

echo "done"
