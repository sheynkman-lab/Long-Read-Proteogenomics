# Peptide Analysis

This module prepares a table comparing mass spec MM peptide results using different databases

## Input

1. gene isoname file: map transcript name to gene name 
2. Gencode peptides file: AllPeptides file from mass spec search using Gencode 
3. Pacbio peptides file: Pacbio refined database fasta file 
4. Pacbio six frame translation: file listing all possible peptides that can be detected per gene 

| argument | description | input module |
|----------|-------------|--------------|
| -gmap    | gene isoname file: map transcript name to gene name  | PG_ReferenceTables |
| -gc_pep  | Gencode peptides file: AllPeptides file from mass spec search using Gencode | MetaMorpheus |
| -pb_pep  | Pacbio peptides file: Pacbio refined database fasta file | MetaMorpheus
| -pb_6frm | Pacbio six frame translation: file listing all possible peptides that can be detected per gene in Pacbio Database | PG_6FrameTranslation |


## Output
Output Tables:
------------------------------------------------------------------------------------------
- table comparing pacbio coverage of Gencode peptide results from MM
------------------------------------------------------------------------------------------


## Soure Module(s)
- PG_ReferenceTables
- MetaMorpheus
- PG_6FrameTranslation

## Target Module(s)
- None

## Dependencies: 
- None

## Threads
- None

## Original Source
- None

## Shell
    python ./peptide_analysis.py \
    -gmap ../../results/PG_ReferenceTables/gene_to_isoname.tsv \
    -gc ../../data/AllPeptides_Gencode.psmtsv \
    -pb ../../data/jurkat_orf_refined.fasta \
    -sft ../../data/pacbio_6frm_database_gene_grouped.fasta \
    -odir ../../results/PG_PeptideAnalysis
