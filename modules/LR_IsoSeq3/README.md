# Iso-Seq 3 processing of CCS reads
The PacBio Iso-Seq method produces high-quality, full-length transcripts and can characterize a whole transcriptome with a single SMRT Cell.

## Input
- jurkat.codethon_toy.ccs.bam (CCS reads)
- NEB_primers.fasta (listing of the adapter sequences)
- hg38.fa (can be retrieved from here: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz)

## Output
- jurkat.collapsed.fasta
- jurkat.collapsed.abundance.txt
Note - there are many other outputs of Iso-Seq 3.

## Soure Module(s)
- (Optional) CCS

## Target Module(s)
- SQANTI3

## Dependencies: 
- None

## Threads
- 40+ cpus

## Original Source
- None

## Shell
```
# create an index
pbindex jurkat.codethon_toy.ccs.bam

module load isoseqenv
lima --isoseq --dump-clips --peek-guess -j 40 jurkat.ccs.bam NEB_primers.fasta jurkat.demult.bam
isoseq3 refine --require-polya jurkat.demult.NEB_5p--NEB_3p.subreadset.xml NEB_primers.fasta jurkat.flnc.bam

# clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
isoseq3 cluster jurkat.flnc.bam jurkat.polished.bam --verbose --use-qvs

# align reads to the genome, takes few minutes (40 core machine)
pbmm2 align hg38.fa jurkat.polished.transcriptset.xml jurkat.aligned.bam --preset ISOSEQ --sort -j 40 --log-level INFO

# collapse redundant reads
isoseq3 collapse jurkat.aligned.bam jurkat.collapsed.gff
```
