# Refined Database Generation
*Description of the module*

## Input

| argument | description | input module |
|----------|-------------|-----------------------|
| --orfs | ORF coordinate input file location | ORFCalling
| --pb_fasta | PacBio fasta sequence input file location | SQUANTI


## Output
| argument | description | 
|----------|-------------|
| --redundant | Output redundant accession file location |
| --combined_tsv | Output combined tsv file location |
| --combined_fasta | Output combined fasta file location |
| --agg_tsv | Output aggregated tsv file location |
| --agg_fasta | Output aggregated fasta file location |

## Soure Module(s)
- None

## Target Module(s)
- None

## Dependencies: 
- None

## Threads
- None

## Original Source
- None

## Shell
    python refine_orf.py \
    --orfs /mnt/shared/ubuntu/session_data/data/output/jurcat_cpat.ORF_called.tsv \
    --pb_fasta /mnt/shared/ubuntu/session_data/data/input/jurkat_corrected.fasta \
    --redundant /mnt/shared/ubuntu/session_data/data/refine/output/jurkat_redundant_accessions.txt \
    --combined_tsv /mnt/shared/ubuntu/session_data/data/refine/output/jurkat_orf_combined.tsv \
    --combined_fasta /mnt/shared/ubuntu/session_data/data/refine/output/jurkat_orf_combined.fasta \
    --agg_tsv /mnt/shared/ubuntu/session_data/data/refine/output/jurkat_orf_aggregated.tsv \
    --agg_fasta /mnt/shared/ubuntu/session_data/data/refine/output/jurkat_orf_aggregated.fasta
