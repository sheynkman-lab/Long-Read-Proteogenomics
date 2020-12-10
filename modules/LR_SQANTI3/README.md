# Long REad SMARTLinkCSS 
*Description of the module*

## Input
- jurkat.collapsed.fasta
- jurkat.collapsed.abundance.txt
- gencode.v35.annotation.gtf
- hg38.fa

## Output
- jurkat.params.txt
- jurkat_classification.txt
- jurkat_corrected.faa
- jurkat_corrected.fasta
- jurkat_corrected.gtf
- jurkat_junctions.txt
- jurkat_sqanti_report.pdf

## Soure Module(s)
- isoseq

## Target Module(s)
- transcriptomeanalysis

## Dependencies: 
- None

## Threads
- None

## Original Source
- None

## Shell
python sqanti3_qc.py jurkat.collapsed.fasta gencode.v35.annotation.gtf hg38.fa --skipORF -o jurkat -d SQANTI3_out/ --fl_count jurkat.collapsed.abundance.txt -n8
