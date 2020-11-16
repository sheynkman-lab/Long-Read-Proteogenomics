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
nextflow run lr_orfcalling.nf --help -profile lr_orfcalling_nextflow.config
```

Then to run a test run on a specific fasta file **ALWAYS USE ABSOLUTE PATHS**
```bash
nextflow run lr_orfcalling.nf --fasta /mnt/shared/ubuntu/session_data/Long-Read-Proteogenomics/data/jurkat_corrected.fasta -profile lr_orfcalling_nextflow.confg
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

### Dependencies: 
- None

### Threads
- None
