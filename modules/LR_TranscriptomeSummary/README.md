# Long Read Transcriptome Summary
*Genomic data is compiled to provide a context for proteomics analysis*

## Input
| argument | description | input module |
|----------|-------------|--------------|
| --sq_out  | classification.txt file from SQUATI Output | SQANTI
| --tpm | Kallisto TPM file location |
| --ribo | Normalized Kallisto Ribodepletion TPM file location |
| --ensg_to_gene | ENSG -> Gene Map file location | PG_ReferenceTables
| --enst_to_isoname | ENST -> Isoname Map file location | PG_ReferenceTables
| --len_stats | Gene Length Statistics table location | PG_ReferenceTables



    parser.add_argument('--odir', '-o', action='store', dest='odir', help='Output Directory')


## Output
| argument | description | output module |
|----------|-------------|--------------|
| --odir | Output Directory |

### Output files
| filename | description | output module |
|----------|-------------|---------------|
| gene_level_tab.tsv | gene level table
| sqanti_isoform_info.tsv | sqanti isoform table |



## Soure Module(s)
- Metamorpheus
- SQANTI3

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
