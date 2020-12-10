# ProteinInferenceAnalysis 
This custom module compares the protein group results from MetaMorpheus searches using different reference databases. The comparison of protein groups can elucidate differences, strengths and weaknesses of the reference databases used.

## Input
- 2 x AllProteinGroups.tsv
- accession mapping file

| argument | description | input module |
|----------|-------------|--------------|
| --pg_fileOne | File location of the AllProteinGroups.tsv file from the first MetaMorpheus search | MetaMorpheus
| --pg_fileTwo | File location of the AllProteinGroups.tsv file from the second MetaMorpheus search | MetaMorpheus |
| --mapping | File location of the Accession mapping file for the converting accessions for comparison purposes | PG_AccessionMapping |


## Output
- excel file containing categories of protein group comparisons (matching, simpler for one model db, partially overlapping or distinct)

## Soure Module(s)
- 2x MetaMorpheus
- AccessionMapping

## Target Module(s)
- None

## Dependencies: 
- Python Packages (pandas, numpy and default dict from collections)

## Threads
- Any

## Original Source
- None

## Shell
- None
