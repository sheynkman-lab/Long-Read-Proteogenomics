
# Deriving the "protein space" (aggregated 6-frame translation) from pacbio data. 


*Description of the module*

## Input
- SQANTI3 isoform structure file
- ENSG to genename map
- PacBio transcript FASTA

## Output
- FASTA of the "protein space" for each gene

## Soure Module(s)
- SQANTI3
- PG_ReferenceTables

## Target Module(s)
- None

## Dependencies: 
- None

## Threads
- None

## Original Source
- None

## Shell
    python make_pacbio6frm_gene_grouped \
    --iso_annot ../SQANTI3_out/jurkat_classification.txt \
    --ensg_gene ../ensg_gene.tsv \
    --sample_fasta ../jurkat_corrected.fasta \
    --output_fasta pacbio_6frm_database_gene_grouped.fasta
