# AccessionMapping
Given several protein databases, map the related isoforms between the databases. Includes the associated genes.

## Input
- Uniprot fasta
- Gencode fasta
- PacBio fasta (from the RefinedDatabase module)

## Output
- an accession mapping file

## Soure Module(s)
- RefinedDatabase

## Target Module(s)
- ProteinInference

## Dependencies: 
- blast
- python package (pandas, numpy, etc.)

## Threads
- 40+ cpus needed for blasting

## Original Source
- None

## Shell
- None
