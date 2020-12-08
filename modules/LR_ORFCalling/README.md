
# Long ORFCalling 
*Description of the module*

This module calls the Open Reading Frames (ORF) with two different ORF calling routines.
TransDecoder and CPAT.

The output of these runs are then filtered by a third program lr_orfcalling.py which cleans up and removes ORFS that are not meaningful.

## Testing locally

To test this part of the pipeline locally install nextflow

To install nextflow - start up a jupyter lab notebook and type
```bash
sudo conda install -c bioconda nextflow
```

To run this pipeline with nextflow - run the test first with displaying the help message
```bash
nextflow run lr_orfcalling.nr --help -profile lr_orfcalling_nextflow.config
```

Then to run a test run on a specific fasta file **ALWAYS USE ABSOLUTE PATHS**
```bash
nextflow run lr_orfcalling.nr --fasta /mnt/shared/ubuntu/session_data/Long-Read-Proteogenomics/data/jurkat_corrected.fasta -profile lr_orfcalling_nextflow.confg
```

## process runTransDecoder

TransDecoder written by Brian Haas, calls the ORF from DNA data obtained from RNA-sequencing
Details may be found here [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki).

It takes two steps - one is to call the longest ORFS with TransDecoder.LongOrfs -t ${params.fasta}.
and second step is predict the best with Transdecoder.Predict -t ${params.fasta}.

All that is required for input is a fasta file.

### Input
- jurkat_corrected.fasta 

### Output
- [outputs explained here](https://github.com/TransDecoder/TransDecoder/wiki)

# ORF calling from CPAT
For each PB transcript and the set of candidate ORFs as determined by CPAT, select the most biologically plausible ORF based on several paramters. These paramters include the coding score and properties of the start ATG.

## Input
- CPAT ORFs (TSV)
- Gencode annotation (GTF)
- SQANTI PB isoform classification
- SQANTI PB annotation (GTF)
- SQANTI PB sequence (FASTA)
- PB to gene map


## Output
- table of best ORF with annotations (TSV)

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
- None

## Script 
    python orf_calling.py \
    --orf_coord ../input/jurkat_cpat.ORF_prob.tsv \
    --gencode ../gencode.v35.annotation.gtf \
    --sample_gtf ../jurkat_corrected.gtf \
    --pb_gene ../pb_to_gene.tsv \
    --classification ../jurkat_classification.txt \
    --sample_fasta ../jurkat_corrected.fasta \
    --output ../output/jurcat_cpat.ORF_called.tsv
