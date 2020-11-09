# ORF calling from CPAT
For each PB transcript and the set of candidate ORFs as determined by CPAT, select the most biologically plausible ORF based on several paramters. These paramters include the coding score and properties of the start ATG.

## Input
- CPAT ORFs (TSV)
- Gencode annotation (GTF)
- SQANTI PB isoform classification
- PB to gene map

## Output
- table of best ORF with annotations

## Soure Module(s)
- CPAT
- SQANTI3

## Target Module(s)
- RefinedDatabase

## Dependencies: 
- pandas

## Threads
- 1-2 cpus

## Original Source
- None

## Shell
- None
