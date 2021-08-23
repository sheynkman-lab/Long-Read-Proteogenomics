# sheynkman-lab/Long-Read-Proteogenomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

NOTE: For an example on where to start counting from in dev phases of development, consult [this informative link](https://www.jering.tech/articles/semantic-versioning-in-practice#semver-in-this-phase)

## v0.1.0

PR: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/pull/103

Initial release of [sheynkman-lab/Long-Read-Proteogenomics](https://github.com/sheynkman-lab/Long-Read-Proteogenomics), created with the [nf-core](https://nf-co.re/) template.

### `Added`
- nf-core based repository structure
- modules and additions to the individual scripts


## v0.2.0

PR: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/pull/132

Enables use of `--config` as a parameter

### `Added`
- Enables use of --config as a parameter
- Enables use of compressed .fa.gz `gencode_translation_fasta` in the process `make_gencode_database`
### `Fixed`
- Fixes test config, add https ZENODO links as env agnostic inputs (runs in cloud, local, cluster, gh-actions)
- Updates initial channel definitions to allow for informative error messages

### `Dependencies`

None added

## (Optional new version)

PR: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/pull/140

Added documentation, scripts to pull data from AWS buckets to Zenodo.

The main reasons why ZENODO vs AWS S3 is a good idea:

1. `Data` versioning (number 1 important reason).  In S3 data can be overwritten for the same path at any point possibly breaking the pipeline.
2. `Cost`, removing data from S3. These data are tiny but the principle stays: The less storage the better
3. `Access`, Academic users familiarity/accessibility - Most reviewers, readers of the pipeline and paper will know ZENODO and will be able to use the data, AWS has an entry barrier for this group of people

### `Added`
- `AWStoZenodo.md` - detail steps for downloading data files from Amazon S3 buckets and upload procedure to `Zenodo`
- `download_aws.sh` - documentation of the files that were extracted from `Amazon S3` buckets.
- `upload_to_zenodo.sh - documentation of the files that were deposited into `Zenodo` repository for use by the `Long-Read-Protegenomics` pipeline.
- `download_and_reconstruct.sh` - script to download files from `Zenodo` for use in execution of the pipeline.

### `Fixed`

None:  To Do:  Modify [`cloudos_jurkat_merged_bam.config`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/conf/cloudos_jurkat_merged_bam.config) file to use the files from Zenodo in the execution of the `main.nf` pipeline

### `Dependencies`

## v0.2.1

PR: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/pull/145

### `Added`
- Optional input into pipeline - metamorpheus_toml for user controlled toml file used in search

### `Fixed`
- novel peptide: set had 'append' rather than 'add' command - bug from change from list to set after incorporating UniProt 

### `Dependencies`

None Added

### `Dependencies`

## v0.2.2

PR: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/pull/153

### `Added`
- output from the full test with running sqanti and running with the two fastq files now on Zenodo with a new version, [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5234651.svg)](https://doi.org/10.5281/zenodo.5234651).  A smaller version of star limited to only chr22 was added to the test repository - all tests run from laptop as well as with GitHub actions.

### `Fixed`
-- all conf files now updated with the new latest version number, test.config (which is redundant with test_with_sqanti.config), test_with_sqanti.config and test_without_sqanti.config'

### `Dependencies`

None added


