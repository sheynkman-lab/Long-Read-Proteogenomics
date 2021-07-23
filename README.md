[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5076056.svg)](https://doi.org/10.5281/zenodo.5076056)
[![reviewdog misspell](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/catch_typos.yml/badge.svg)](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/catch_typos.yml)
[![Testing for Long Reads Proteogenomics](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/ci.yml/badge.svg)](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/ci.yml)
# Sheynkman-Lab/Long-Read-Proteogenomics

`Updated: 2021 July 11`

This is the repository for the **Long-Read Proteogenomics** workflow.  Written in [`Nextflow`](https://www.nextflow.io/), it is a modular workflow beneficial to both the `Transcriptomics` and `Proteomics` fields. The data from both `Long-Read IsoSeq sequencing` with `PacBio` and `Mass spectrometry-based proteomics` used in the classification and analysis of protein isoforms expressed in `Jurkat` cells and described in the publication `Enhanced protein isoform characterization through long-read proteogenomics`, which will be made public in Fall 2021.

A goal in the biomedical field is to delineate the protein isoforms that are expressed and have pathophysiological relevance. Towards this end, new approaches are needed to detect protein isoforms in clinical samples. Mass spectrometry (MS) is the main methodology for protein detection; however, poor coverage and incompleteness of protein databases limit its utility for isoform-resolved analysis. Fortunately, long-read RNA-seq approaches from PacBio and Oxford Nanopore platforms offer opportunities to leverage full-length transcript data for proteomics.

We introduce enhanced protein isoform detection through integrative “long read proteogenomics”. The core idea is to leverage long-read RNA-seq to generate a sample-specific database of full-length protein isoforms. We show that incorporation of long read data directly in the MS protein inference algorithms enables detection of hundreds of protein isoforms intractable to traditional MS. We also discover novel peptides that confirm translation of transcripts with retained introns and novel exons. Our pipeline is available as an open-source Nextflow pipeline, and every component of the work is publicly available and immediately extendable.

Proteogenomics is providing new insights into cancer and other diseases. The proteogenomics field will continue to grow, and, paired with increases in long-read sequencing adoption, we envision use of customized proteomics workflows tailored to individual patients.

We acknowledge the beginning kernels of this work were formed during the Fall of 2020 at the [`Cold Spring Harbor Laboratory Biological Data Science Codeathon`](https://datascience.nih.gov/news/cold-spring-harbor-laboratory-biological-data-science-codeathon).  

We acknowledge Lifebit and the use of their platform Lifebit's CloudOS key in development of the open source software Nextflow workflow used in this work.  
<p align="center"><img src="https://github.com/lifebit-ai/dry-bench-skills-for-researchers/blob/adds-mini-courses/assets/lifebitCloudOS.png"  width="250" align="right" ></p>

## How to use this repository and Quick Start

This workflow was used end-to-end in with the Nextflow main.nf workflow.  However, processes within this workflow maybe useful for others and are self-contained with clear identification of inputs and outputs.  These are documented on the [`wiki`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki) within this repository.

You can clone this GitHub repository and run this on an appropriately sized machine.  

## Quick Start

### `1. Create the environment`

Always best practice to work within a controlled environment.  We use [`conda`] to create this environment.

```bash
conda init bash
exec -l bash
```
Create and activate a new conda environment `lrp`.

```bash
conda create -n lrp
conda activate lrp
```
### `2. Download nextflow`

Download Nextflow (use [`anaconda search`] for the correct syntax.  Searching shows [`Nextflow`] is available on the `Bioconda` channel.

```bash
conda install -c bioconda nextflow -y
```

### `3. Clone this repository`

Now with the environment ready - you can either clone the environment or run the repository without cloning. See [`Wiki`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki) for all of the available options when running the pipeline.

Optionally, clone the [Long-Read-Proteomics](https://github.com/sheynkman-lab/Long-Read-Proteogenomics) repository.

```bash
git clone https://github.com/sheynkman-lab/Long-Read-Proteogenomics.git
```

### `4. Execute the workflow as a test`
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5109695.svg)](https://doi.org/10.5281/zenodo.5109695)

One test that can be run uses a test configuration file stored on [`Zenodo`](https://doi.org/10.5281/zenodo.5109695).

Running a test from the clone can be done with the command - without metamorphesis

```bash
nextflow run main.nf -profile test
```
This was done on a MacBook Pro running BigSur Version 11.4 with 16 GB 2667 MHz DDR48 RAM and a 2.3 GHz 8-Core Intel Core i9 processor.  

## Documentation and Pipeline Vignette

The sheynkman-lab/Long-Read-Proteogenomics pipeline comes with details about each of the processes that make up the pipeline are found in the [`Wiki`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki).   In this you will find:

1. [`Third-party tools`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/3rd-Party-Pipeline-Tools)
2. [`Input parameters`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Input-Parameters)
3. [`Output files`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Output-Files)
4. [`Pipeline processes descriptions`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Pipeline-Processes)
5. [`Pipeline vignette`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Pipeline-Vignette)

## Pipeline overview

The pipeline accepts as input raw PacBio data and performs the assembly of predicted protein isoforms with high probability of existing in the sample. This database is then used in [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to search raw mass spectrometry data against the PacBio reference. MetaMorpheus will use protein isoform read counts during protein inference. Two other protein databases are employed for the purposes of comparison. One is from [UniProt](https://www.uniprot.org/) and the other is from [GENCODE](https://www.gencodegenes.org/). A series of [Jupyter notebooks](https://github.com/sheynkman-lab/LRPG-Manuscript) can be used to perform all final comparisons and data analysis. 

![LRP Pipeline_v2](https://user-images.githubusercontent.com/12956799/110397183-5c041b00-803f-11eb-9ba7-02352dab5656.png)

## Using Zenodo

To make the data more accessible and FAIR, the indexed files were transferred to Zenodo using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload) from the `University of Virginia's Gloria Sheynkman Lab` Amazon `S3` buckets.

The data were prepared and stored in the development of the `Sheynkman Lab Long-read Proteogenomics Pipeline`.

Using Nextflow, configuration items can access locations in Google Compute Platform (GCP) buckets (`gs://`), Amazon Web Services (AWS) buckets (`s3://`) and Zenodo locations (`https://`) seamlessly.

The main reasons why ZENODO vs AWS S3: or GCP GS: are:

1. `Data versioning` (of primary importance): In S3 or GS buckets, data can be overwritten for the same path at any point, possibly breaking the pipeline.
2. `Cost`: These datasets are tiny but the principle stays: The less storage the better
3. `Access`: Most users of the pipeline can most easily access `ZENODO` and will be able to use the data. AWS and GCP has an entry barriers.

Details on how these data were transferred and moved from `AWS S3:` buckets are described in the [`AWS to Zenodo`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/AWStoZenodo.md) markdown document in this repository.

## Contributors

- [Christina Chatzipantsiou](https://orcid.org/0000-0002-4257-7241)
- [Benjamin Jordan](https://orcid.org/0000-0003-2268-5226)
- [Simran Kaur](https://orcid.org/0000-0003-1536-5996)
- [Raymond Leclair](https://orcid.org/0000-0002-7233-6588)
- [Anne Deslattes Mays](https://orcid.org/0000-0001-7951-3439)
- [Madison Mehlferber](https://orcid.org/0000-0001-9598-8100)
- [Rachel M. Miller](http://orcid.org/0000-0003-1461-6386)
- [Robert J. Millikin](https://orcid.org/0000-0001-7440-3695)
- [Kyndalanne Pike](https://orcid.org/0000-0002-5906-2340)
- [Gloria M. Sheynkman](https://orcid.org/0000-0002-4223-9947)
- [Michael R. Shortreed](https://orcid.org/0000-0003-4626-0863)
- [Isabella Whitworth](https://orcid.org/0000-0001-9753-5666)

This is a joint project between the [Sheynkman Lab](https://med.virginia.edu/sheynkman-lab/), the [Smith Lab](https://smith.chem.wisc.edu/), [Lifebit](https://lifebit.ai/) and Science and Technology Consulting, LLC.


## Repository template

This pipeline was generated using a modification of the nf-core template. 
You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
