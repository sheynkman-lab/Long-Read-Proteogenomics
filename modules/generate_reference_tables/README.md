# Reference Table Generation
*Description of the module*

## Input
| command | description |
|---------| ----------- |
| --gtf | Gencode GTF input file location
| --fa | Gencode Fafsa input file location

## Output

| command | description |
|--------| ---------- |
| --ensg_gene | ensg to gene output file location |
|--enst_isoname | enst to isoname output file location |
|--gene_ensp | Gene to ensp output file location |
|--gene_isoname | Gene to isoname output file location |
|--isoname_lens | Isoname length table output location |
|--gen_lens | Gene Length statistics output location |

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
python prepare_reference_tables \
--gtf ../input/gencode.gtf \
--fa ../input/gencode.fafsa \
--ensg_gene ../results/ensg_gene.tsv \
--enst_isoname ../results/enst_isoname.tsv \
--gene_ensp ../reesults/gene_ensp.tsv \
--gene_isoname ../results/gene_isoname.tsv \
--isoname_lens ../results/isoname_lens.tsv \
--gene_lens ../results/gene_lens.tsv