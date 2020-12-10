
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