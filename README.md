# sheynkman-lab/Long-Read-Proteogenomics

**A workflow for Delineating the Human Proteome at Isoform Resolution by Integration of Long-read Proteogenomics and Mass Spectrometry**.

[Cold Spring Harbor Laboratory Biological Data Science Codeathon](https://datascience.nih.gov/news/cold-spring-harbor-laboratory-biological-data-science-codeathon)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4256568.svg)](https://doi.org/10.5281/zenodo.4256568)

## Quick Start

1. Make sure you have [`nextflow`](https://nf-co.re/usage/installation), [`Docker`](https://docs.docker.com/engine/installation/) installed
3. Start running your own analysis!

    ```bash
    nextflow run sheynkman-lab/Long-Read-Proteogenomics -profile test
    ```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The sheynkman-lab/Long-Read-Proteogenomics pipeline comes with documentation about the pipeline: [usage](docs/usage.md) and [output](docs/output.md).

### Pipeline overview

The pipeline accepts as input raw Pac-Bio data and performs the assembly of an accurate list of protein isoforms with high probability of existing in the sample. This database is then used in [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus) to search raw mass spectrometry data against the Pac-Bio reference. MetaMorpheus will use protein isoform read counts during protein inference. Two other protein databases are employed for the purposes of comparison. One is from [UniProt](https://www.uniprot.org/) and the other is from [GENCODE](https://www.gencodegenes.org/). A Jupyter notebook performs all final comparisons and data analysis. 

![pipeline](https://user-images.githubusercontent.com/16841846/98399434-fa4b8a00-2027-11eb-953b-edb440c7ff8e.png)


## Background

Protein isoforms are the direct translational product of fully spliced mRNA molecules. Protein isoforms can be modified during or subsequent to translation with additional chemical moities (e.g. phosphorylation or acetylation) or they can be cleaved resulting in a [proteoform](https://www.nature.com/articles/nmeth.2369), which is the ultimate biological actor in many important biological processes. At a high level, protein isoforms can be predicted from genomic sequencing data and then observed by mass spectrometry. Despite impressive technological achievements in both realms (sequencing and mass spectrometry), many gaps exist in our ability to comprehensively identify all protein isoforms even for a single sample. Scientists ability to accomplish this goal depends on having detailed an accurate knowledge of all protein coding mRNA isoforms, comprehensive mass spectrometry data covering at least one unique region of each protein isoform, and a protein inference algorithm that can faithfully and accurately attribute observed peptides to the proper parent isoform. We provide below an overiew of the key remaining challenges and then provide for the first time a pipeline that solves these challenges.

### Protein Isoforms

Knowledge of a full-length transcriptome can provide for an empirically-derived predicted set of protein isoforms, which can serve as accurate and more precise models for protein inference. Third generation sequencing, exemplified recently by [Pac-Bio](https://www.pacb.com/) can, for the first time, shed light on full-length protein isoforms. Until now, protein isoforms were inferred through transript reconstruction on next generation sequencing data. However, this was a frought process with many errors. With the advent of long-read sequencing, we can observe full-length, fully-spliced mRNA transcripts that can be translated into protein sequencing for use in subsequent mass spectrometry experiments. A major remaining challenge is the identification of all open reading frames (ORFs). 

### Mass Spectrometry Data

Bottom-up mass spectrometry is the premier method for protein identification. Mass-spectrometry, as as technology, provides a means to rapidly identify peptides produces by proteolytic digestion of intact proteins isoforms. It is fast and sensitive. Well done experiments frequently identify as many as 10,000 proteins in a single analysis. Yet, much can be done to improve the depth and accuracy of such experiments, especially for comprehensive identification of protein isoforms. First and foremost, the dominant choice of protease for bottom-up mass spectrometry is trypsin. Trypsin digest whole proteins into manageable peptides that are easily separated by HPLC and identified by mass spectrometry. However, identification of a protein isoform requires at minimum a single peptide that can be uniquely ascribed to that isoform. Here, trypsin alone simply cannot deliver enough unique peptides to identify all protein isoforms in a sample. The reason is that many tryptic peptides are too short or too long for effective mass-spec analysis. In addition, many tryptic peptides are shared between many protein isoforms giving them litte informative value. Recently, [Miller](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00330) demonstrated that use of multiple proteases for a single sample, can greatly improve protein inference by significantly increasing the number of unique peptides detected. Frequently, protein isoforms can have multiple unique peptides for added identification confidence. 

### Protein Isoform Inference

[Protein inference](https://www.sciencedirect.com/science/article/abs/pii/S187439191630344X?via%3Dihub) is the process of guessing which proteins are present in a sample based on limited peptide evidence. Bottom-up proteomics, by definition, deals only in peptides, which are the pieces of a protein available for analysis following digestion with a protease. Top-down proteomics would be the preferred method for protein isoform detections because it analyzes intact proteoforms. However, at the present time, it lacks the sensitivity that bottom-up has, yielding only fractional proteome coverage. In bottom-up, a key challenge is taking all of the identified peptides and then attempting to use them to infer presence of the protein isoforms from which they were derived. This process is aided greatly by deeper coverage of peptides unique to each isoform in the sample. Still it is not a solved problem. Here, in this project, we will integrate protein isoform presence as measured by copy number from the Pac-Bio data as a Bayesian prior in the protein inference algorithm. 

## Contributors

- [Christina Chatzipantsiou](https://orcid.org/0000-0002-4257-7241)
- [Benjamin Jordan](https://orcid.org/0000-0003-2268-5226)
- [Simran Kaur](https://orcid.org/0000-0003-1536-5996)
- [Raymond Leclair](https://orcid.org/0000-0002-7233-6588)
- [Anne Deslattes Mays](https://orcid.org/0000-0001-7951-3439)
- Madison Mehlferber
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
