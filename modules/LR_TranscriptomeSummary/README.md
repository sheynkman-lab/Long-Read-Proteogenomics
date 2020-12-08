# Long Read Transcriptome Summary
*Genomic data is compiled to provide a context for proteomics analysis*

## Input
| argument | description |
|----------|-------------|
| gtf_file     | gencode gtf file 
| squanti_out  | classification.txt file from SQUATI Output 
| fa_file      | gencode fafsa file
| tpm_file     |
| yangpolyA    | file prepared manually
| ribodep_tpm           | expects normalized data
| pbacc_to_gene_file    | mapping file from Gloria
| mm_out                | Metamorpheus outputs 

## Output
- squanti isoform table
- gene level table
- gene level and proteomic comparison table 

## Soure Module(s)
- Metamorpheus
- SQUANTI3

## Target Module(s)
- None

## Dependencies: 
- numpy
- pandas
- pathlib from Path

## Threads
- None

## Original Source
- None

## Shell
    python transcriptome_summary.py \
    --sq_out ../input/sq_out.tsv \
    --tpm ../input/sq_out.tsv \
    --ribo ../input/sq_out.tsv \
    --ensg_to_gene ../input/sq_out.tsv \
    --enst_to_isoname ../input/sq_out.tsv \
    --len_stats ../input/sq_out.tsv \
    --odir ../output
