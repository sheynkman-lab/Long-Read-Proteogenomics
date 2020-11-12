# Long REad SMARTLinkCSS 
*Description of the module*

To test this part of the pipeline locally install nextflow

To install nextflow - start up a jupyter lab notebook and type
```bash
sudo conda install -c bioconda nextflow
```

To run this pipeline with nextflow - run the test first with displaying the help message
```bash
nextflow run lr_orfcalling.nr --help -profile lr_orfcalling_nextflow.config
```

Then to run a test run on a specific fasta file ALWAYS USE ABSOLUTE PATHS
```bash
nextflow run lr_orfcalling.nr --fasta /mnt/shared/ubuntu/session_data/Long-Read-Proteogenomics/obs
clear

data/jurkat_corrected.fasta -profile lr_orfcalling_nextflow.confg
```

## process runTransDecoder

TransDecoder written by Brian Haas, calls the ORF from DNA data obtained from RNA-sequencing
Details may be found here [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki).

It takes two steps - one is to call the longest ORFS with TransDecoder.LongOrfs -t ${params.fasta}.
and second step is predict the best with Transdecoder.Predict -t ${params.fasta}.

All that is required for input is a fasta file.

## Input
- jurkat_corrected.fasta 

## Output
- [outputs explained here](https://github.com/TransDecoder/TransDecoder/wiki)
- D

## Dependencies: 
- None

## Threads
- None

