# Run CPAT on PacBio transcript file
CPAT finds the most likely ORF from a full-length transcript, and outputs all possible ORFs with accompanying metrics.

## Input
- Human_Hexamer.tsv
- Human_logitModel.RData
- jurkat_corrected.fasta (from SQANTI3 module)

## Output
- Table of orf calls and their scores (e.g., jurkat_cpat.ORF_prob.tsv)

## Soure Module(s)
- SQANTI3

## Target Module(s)
- LR_ORFCalling
- RefinedDatabase

## Dependencies: 
- CPAT

## Threads
- 1-2 cores

## Original Source
- https://cpat.readthedocs.io/en/latest/

## Shell
- None
