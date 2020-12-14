
# Long Read ORFCalling 
*Description of the module*

This module calls the Open Reading Frames (ORF) based on candidates generated from CPAT.

The output of these runs are then filtered by a program orfcalling.py which cleans up and removes ORFS that are least biologically plausible.

## Testing locally

To test this part of the pipeline locally install nextflow

To install nextflow - start up a jupyter lab notebook and type
```bash
sudo conda install -c bioconda nextflow
```

To run this pipeline with nextflow - run the test first with displaying the help message
```bash
nextflow run lr_orfcalling.nf --help -profile nextflow.config
```

Then to run a test run on a specific fasta file **ALWAYS USE ABSOLUTE PATHS**
```bash
nextflow run lr_orfcalling.nf --fasta /mnt/shared/ubuntu/session_data/Long-Read-Proteogenomics/data/jurkat_corrected.fasta -profile lr_orfcalling_nextflow.confg
```

# ORF calling from CPAT
For each PB transcript and the set of candidate ORFs as determined by CPAT, select the most biologically plausible ORF based on several paramters. These paramters include the coding score and properties of the start ATG.

## Input
| argument | description | input module |
|----------|-------------|--------------|
| --orf_coord |CPAT ORFs (TSV) | LR_CPAT |
| --gencode | Gencode annotation (GTF) | 
| --sample_gtf | SQANTI PB annotation (GTF) | SQANTI |
| --pb_gene | PacBio : Genecode Cross Ref | PG_ReferenceTables
| --classification | SQANTI PB isoform classification | SQANTI |
| --sample_fasta | SQANTI PB sequence (FASTA) | SQANTI |

## Output
| argument | description | output module |
|----------|-------------|--------------|
| --output | table of best ORF with annotations (TSV) | RefineDatabaseGeneration

## Soure Module(s)
- CPAT
- SQANTI3

## Target Module(s)
- RefinedDatabase

## Dependencies: 
- pandas

## Threads
- 1-2 cpus

### Dependencies: 
- None

### Threads
- 1-2 cpus

## Script 
    python orf_calling.py \
    --orf_coord ../input/jurkat_cpat.ORF_prob.tsv \
    --gencode ../gencode.v35.annotation.gtf \
    --sample_gtf ../jurkat_corrected.gtf \
    --pb_gene ../pb_to_gene.tsv \
    --classification ../jurkat_classification.txt \
    --sample_fasta ../jurkat_corrected.fasta \
    --output ../output/jurkat_cpat.ORF_called.tsv
