# SQANTI3 analysis
Isoform annotations

## Input
- jurkat.collapsed.gff
- jurkat.collapsed.abundance.txt
- gencode.v35.annotation.gtf
- hg38.fa (must be uncompressed)

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

```bash
sqanti3_qc.py \
jurkat.collapsed.canonical.chr22.gff \
gencode.v35.annotation.chr22.gtf \
GRCh38.primary_assembly.genome.fa \
--skipORF \
-o jurkat \
-d SQANTI3_out/ \
--fl_count jurkat.collapsed.abundance.txt  \
--gtf
```

### NOTES:

#### Using the docker image

When running the bash command we can use the docker image for acquiring the necessary dependencies. To mount the container and use the input files in the current working directory we can run the following command:

```
docker run --rm --entrypoint /bin/bash  -v $PWD:$PWD -w $PWD -it  gsheynkmanlab/sqanti3
conda activate sqanti3
```

#### Non canonical chromosomes

If there are entries of non canonical chromosomes in the eg. `jurkat.collapsed.chr22.gff` input file, make sure you are using the required reference file or preprocess the gff file by removing `chr22_*` lines.

to do so:

```bash
grep -v chr22_ jurkat.collapsed.chr22.gff > jurkat.collapsed.canonical.chr22.gff
```

#### Compressed fasta file


<details>
<summary>
The following error was encountered when using the compressed fasta file:
</summary>


```console
Reading genome fasta /home/ec2-user/Long-Read-Proteogenomics/modules/LR_SQANTI3/GRCh38.primary_assembly.genome.fa.gz....
Traceback (most recent call last):
  File "/opt/sqanti3/sqanti3_qc.py", line 2291, in <module>
    main()
  File "/opt/sqanti3/sqanti3_qc.py", line 2275, in main
    run(args)
  File "/opt/sqanti3/sqanti3_qc.py", line 1701, in run
    genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.genome), 'fasta'))
  File "/opt/sqanti3/sqanti3_qc.py", line 1701, in <genexpr>
    genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.genome), 'fasta'))
  File "/opt/conda/envs/sqanti3/lib/python3.7/site-packages/Bio/SeqIO/FastaIO.py", line 186, in FastaIterator
    for title, sequence in SimpleFastaParser(handle):
  File "/opt/conda/envs/sqanti3/lib/python3.7/site-packages/Bio/SeqIO/FastaIO.py", line 48, in SimpleFastaParser
    for line in handle:
  File "/opt/conda/envs/sqanti3/lib/python3.7/codecs.py", line 322, in decode
    (result, consumed) = self._buffer_decode(data, self.errors, final)
```
    
</details>

Make sure you are using the uncompressed version.


#### `-n8` (splits)

Omitted for now, as this returns an error message, as follows

