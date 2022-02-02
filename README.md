[![reviewdog misspell](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/catch_typos.yml/badge.svg)](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/catch_typos.yml)[![Testing for Long Reads Proteogenomics without Sqanti](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/ci_test_without_sqanti.yml/badge.svg)](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/ci_test_without_sqanti.yml)

This Repository contains the complete software and documentation to execute the Long-Read-Proteogenomics Workflow.

## Digital Object Identifiers

For the Genome Biology Manuscript.
| DOI    | Description       |
| ------------- | --------------------------------------------------------------------------- |
| <a href="https://doi.org/10.5281/zenodo.5920817" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5920817.svg?raw=true" alt="drawing" width="500"/></a> | Contains the version of the repository used for execution and generation of data |
| <a href="https://doi.org/10.5281/zenodo.5703754" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5703754.svg?raw=true" alt="drawing" width="500"/></a> | Contains the input data from Jurkat Samples and Reference data used in execution of the Long-Read-Proteogenomics workflow contained in this repository| 
| <a href="https://doi.org/10.5281/zenodo.5920920" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5920920.svg?raw=true" alt="drawing" width="500"/></a> | Contains the output data from executing the Long-Read-Proteogenomics workflow using the Zenodo version of this repository |
|  <a href="https://doi.org/10.5281/zenodo.5920847" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5920847.svg?raw=true" alt="drawing" width="500"/></a> | Contains the version of analysis codes and codes for generating the figures using as input the output data from executing the Long-Read-Proteogenomics Workflow version specified above |
|  <a href="https://doi.org/10.5281/zenodo.5234651" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5234651.svg?raw=true" alt="drawing" width="500"/></a> | Contains the Test Data used with the GitHub Actions to ensure changes to this repository still execute and perform correctly |

| Sequence Read Archive (SRA) Project Reference    | Description       |
| ------------- | --------------------------------------------------------------------------- |
| [PRJNA783347](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA783347) | Long-Read RNA Sequencing Project for Jurkat Samples |
| [PRJNA193719](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA193719) | Short-Read RNA Sequencing Project for Jurkat Samples |


# Sheynkman-Lab/Long-Read-Proteogenomics

`Updated: 2022 January 30`

This is the repository for the **Long-Read Proteogenomics workflow**.  Written in [`Nextflow`](https://www.nextflow.io/), it is a modular workflow beneficial to both the `Transcriptomics` and `Proteomics` fields. The data from both `Long-Read IsoSeq sequencing` with `PacBio` and `Mass spectrometry-based proteomics` used in the classification and analysis of protein isoforms expressed in `Jurkat` cells and described in the publication `Enhanced protein isoform characterization through long-read proteogenomics`, which will be made public in Fall 2022.

The output data resulting from the execution of this workflow for the [Manuscript](). May be found here [insert Zenodo Reference here].   The Analysis to produce the figures for the manuscript may be found in the companion repository [**Long-Read Proteogenomics Analysis**](https://github.com/sheynkman-lab/Long-Read-Proteogenomics-Analysis)

A goal in the biomedical field is to delineate the protein isoforms that are expressed and have pathophysiological relevance. Towards this end, new approaches are needed to detect protein isoforms in clinical samples. Mass spectrometry (MS) is the main methodology for protein detection; however, poor coverage and incompleteness of protein databases limit its utility for isoform-resolved analysis. Fortunately, long-read RNA-seq approaches from PacBio and Oxford Nanopore platforms offer opportunities to leverage full-length transcript data for proteomics.

We introduce enhanced protein isoform detection through integrative “long read proteogenomics”. The core idea is to leverage long-read RNA-seq to generate a sample-specific database of full-length protein isoforms. We show that incorporation of long read data directly in the MS protein inference algorithms enables detection of hundreds of protein isoforms intractable to traditional MS. We also discover novel peptides that confirm translation of transcripts with retained introns and novel exons. Our pipeline is available as an open-source Nextflow pipeline, and every component of the work is publicly available and immediately extendable.

Proteogenomics is providing new insights into cancer and other diseases. The proteogenomics field will continue to grow, and, paired with increases in long-read sequencing adoption, we envision use of customized proteomics workflows tailored to individual patients.

We acknowledge the beginning kernels of this work were formed during the Fall of 2020 at the [`Cold Spring Harbor Laboratory Biological Data Science Codeathon`](https://datascience.nih.gov/news/cold-spring-harbor-laboratory-biological-data-science-codeathon).  

We acknowledge Lifebit and the use of their platform Lifebit's CloudOS key in development of the open source software Nextflow workflow used in this work.  
<p align="center"><img src="https://github.com/lifebit-ai/dry-bench-skills-for-researchers/blob/main/assets/lifebitCloudOS.png"  width="250" align="right" ></p>

## How to use this repository and Quick Start

This workflow is complex, bringing together two measurement technologies in a long-read proteogenomics approach for integrating sample-matched long-read RNA-seq and MS-based proteomics data to enhance isoform characterization.   To orient the user with the steps involved in the transformation of raw measurement data to these fully resolved, identified and annotated results, we have developed this quick start, wiki documentation including vignettes.  

### How to use this repository

This repository is organized into modules and parts of this repository could be useful to different researchers to annotate their own raw data.   The workflow is written in [`Nextflow`](https://www.nextflow.io/), allowing it to be run on virtually any platform with alterations to the configurations and other adaptations.   The visitor is encourated to fork clone and adapt and contribute.   All are encouraged to use [`GitHub Issues`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/issues) to communicate with the contributors to this open source software project.   Software addtions, modifications and contributions are done through [`GitHub Pull Requests`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/pulls)  

Module processes details are documented within the [`Wiki`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki) within this repository.  As well as linked to the third party resources used in this workflow.

Vignettes have been developed to go into greater detail and walk the visitor through the [`visualization capabilities of the final annotated results`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Vignette---How-to-Visualize-Results-from-Jurkat-Analysis) and to walk the visitor through the workflow with presented here with the [`quick start`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Vignette-Long-Read-Proteogenomics-Workflow-with-Test-Data)

### Quick Start

This quick start and steps were performed on a MacBook Pro running BigSur Version 11.4 with 16 GB 2667 MHz DDR48 RAM and a 2.3 GHz 8-Core Intel Core i9 processor.

The visitor will be walked through the pre-requisites, clone the library and execute with demonstration data also used in the [`GitHub Actions`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/actions/workflows/ci_test_without_sqanti.yml). 

#### Obtain the Desktop DockerHub Application

<p align="center"><img src="https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/docs/images/Moby-logo.png"  width="100" align="right" ></p>

In this quick start, [`Dockerhub Desktop Application for the Mac with an Intel Chip`](https://hub.docker.com/editions/community/docker-ce-desktop-mac) was used.
Follow the instructions there to install.

#### Configure the Desktop DockerHub Application

On the MacBook Pro running BigSur Version 11.4 with 16 GB Ram, It was necessary to configure the Dockerhub resources to use **`6GB`** of Ram.

<img src="https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/docs/images/DockerHubDesktopResourceConfigurationMacWithIntelChip.png"  width="600" height="400">

#### Obtain and install miniconda

On the MacBook Pro, the [`64-bit version of miniconda was downloaded and installed`](https://docs.conda.io/en/latest/miniconda.html) follow the installation instructions.

#### Create and activate a new conda environment `lrp`.

To begin, open a terminal window, ensuring the miniconda installation has completed, reboot the terminal shell.
On the Mac, this is done within a `zsh` shell environment.

```bash
exec -l zsh
```

If you already have the environment, you can see what conda environments you have with the following commnad:
```bash
conda info --envs
```

If you haven't already created a conda environment for this work, create and activate it now.
```bash
conda create -n lrp
conda activate lrp
```

#### Install Nextflow.

Install and set the Nextflow version.

```bash
conda install -c bioconda nextflow -y
export NXF_VER=20.01.0
```

#### Clone this repository 

Now with the environment ready, we can clone.

```bash
git clone https://
.com/sheynkman-lab/Long-Read-Proteogenomics
cd Long-Read-Proteogenomics
```

#### Run the pipeline with the test_without_sqanti.config

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5234651.svg)](https://doi.org/10.5281/zenodo.5234651)

This Quick start uses the `test_without_sqanti.config` configuration file found in the `conf` directory of this repository.

```bash
nextflow run main.nf --config conf/test_without_sqanti.config 
```

For details regarding the processes and results produced, please see the [`Wiki`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki) and the  [`Vignette: Workflow with test data`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Vignette-Long-Read-Proteogenomics-Workflow-with-Test-Data).

To visualize results, please see the [`visualization capabilities of the final annotated results`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Vignette---How-to-Visualize-Results-from-Jurkat-Analysis).

## Documentation and Workflow Vignettes

The sheynkman-lab/Long-Read-Proteogenomics pipeline comes with details about each of the processes that make up the pipeline are found in the [`Wiki`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki).   In this you will find:

1. [`Third-party tools`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/3rd-Party-Pipeline-Tools)
2. [`Input parameters`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Input-Parameters)
3. [`Output files`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Output-Files)
4. [`Pipeline processes descriptions`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Pipeline-Processes)
5. [`Vignette: Visualization`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Vignette---How-to-Visualize-Results-from-Jurkat-Analysis)
6. [`Vignette: Workflow with test data`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/wiki/Vignette-Long-Read-Proteogenomics-Workflow-with-Test-Data)

## Workflow overview

The workflow accepts as input raw PacBio data and performs the assembly of predicted protein isoforms with high probability of existing in the sample. This database is then used in [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to search raw mass spectrometry data against the PacBio reference. MetaMorpheus will use protein isoform read counts during protein inference. Two other protein databases are employed for the purposes of comparison. One is from [UniProt](https://www.uniprot.org/) and the other is from [GENCODE](https://www.gencodegenes.org/). A series of [Jupyter notebooks](https://github.com/sheynkman-lab/LRPG-Manuscript) can be used to perform all final comparisons and data analysis. 

![LRP Pipeline_v2](https://user-images.githubusercontent.com/12956799/110397183-5c041b00-803f-11eb-9ba7-02352dab5656.png)

## Using Zenodo

To make the data more accessible and FAIR, the indexed files were transferred to Zenodo using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload) from the `University of Virginia's Gloria Sheynkman Lab` Amazon `S3` buckets. 

Using Nextflow, configuration items can access locations in Google Compute Platform (GCP) buckets (`gs://`), Amazon Web Services (AWS) buckets (`s3://`) and Zenodo locations (`https://`) seamlessly.

The main reasons why ZENODO vs AWS S3: or GCP GS: are:

1. `Data versioning` (of primary importance): In S3 or GS buckets, data can be overwritten for the same path at any point, possibly breaking the pipeline.
2. `Cost`: These datasets are tiny but the principle stays: The less storage the better
3. `Access`: Most users of the pipeline can most easily access `ZENODO` and will be able to use the data. AWS and GCP has an entry barriers.

Details on how these data were transferred and moved from `AWS S3:` buckets are described in the [`AWS to Zenodo`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/AWStoZenodo.md).

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
