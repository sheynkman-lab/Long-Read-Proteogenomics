[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5076056.svg)](https://doi.org/10.5281/zenodo.5076056)

# Using Zenodo

To make the data more accessible and FAIR, the indexed files were transfered to Zenodo using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload) from the `University of Virginia's Gloria Sheynkman Lab` Amazon `S3` buckets.

The data were prepared and stored in the development of the [`Sheynkman Lab Long-read Proteogenomics Pipeline`](https://github.com/sheynkman-lab/Long-Read-Proteogenomics).

Using Nextflow, configuration items can access locations in Google Compute Platform (GCP) buckets (`gs://`), Amazon Web Services (AWS) buckets (`s3://`) and Zenodo locations (`https://`) seamlessly.

The main reasons why ZENODO vs AWS S3: or GCP GS: are:

1. `Data versioning`: (number 1 important reason), In S3 or GS, data can be overwritten for the same path at any point possibly breaking the pipeline.
2. `Cost`: These data are tiny but the principle stays: The less storage the better
3. `Access`: Most reviewers, readers of the pipeline and paper will know `ZENODO` and will be able to use the data, AWS and GCP has an entry barrier for many.

## Prepare your environment

<p align="center"><img src="https://github.com/lifebit-ai/dry-bench-skills-for-researchers/blob/main/assets/lifebitCloudOS.png"  width="250" align="right" ></p>

We used [`JupyterLab`](https://jupyter.org/) terminal windows within the [`Lifebit's CloudOS`](https://lifebit.ai/cloudos) application.  The documentation here is to be transparent about the data that were moved from the Amazon S3 buckets to Zenodo.  Similar steps may be performed on other datasets following the same step-by-step procedure in an appropriate sized terminal, running a unix environment.

#### i. initialize the bash environment

```bash
conda init bash
exec -l bash
```

#### ii. create and activate a new conda environment `lrp`.
```bash
conda create -n lrp
conda activate lrp
```

Now with the environment ready - begin the pre-requisites

## Clone the [zenodo-upload](https://github.com/jhpoelen/zenodo-upload) repository

```bash
git clone https://github.com/jhpoelen/zenodo-upload.git
```

## Satisfy the Pre-requisites for [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload)

Following the instructions from [`README.md`](https://github.com/jhpoelen/zenodo-upload#readme).
The following steps were performed.

#### i. Install [`jq`](https://stedolan.github.io/jq/).  

```bash
conda install -c conda-forge jq -y
```

#### ii. Installed `curl`

```bash
conda install -c conda-forge curl -y
```

#### iii. Bash

already satisfied.

## Preparing to use [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload)

## Transfer data from an S3 bucket

Data were on aws.  Add `awscli` the command line library for AWS S3 buckets.
For organizational principals, begin in the data directory of the recently cloned file.

#### i. install command-line interface to Amazon Web Services `awscli`

```bash
conda install -c conda-forge awscli -y
```

#### ii. clone the repository

```bash
git clone https://github.com/sheynkman-lab/Long-Read-Proteogenomics.git
```

#### iii. download the files

```bash
cd Long-Read-Proteogenomics/data
bash ../download_aws.sh
```

## Prepare data for upload to `Zenodo`
The download script has flattened the hierarchy of files.  Ensuring the scripts run appropriately, we need to recreate the directory structure where many of these files were previously arranged into folders.  These folders will be replaced by tar'd and zipped files that are stored in Zenodo. 

We do this with tar (system provided) and [pigz](http://zlib.net/pigz/)

### tar and pigz

For each previously arranged file that was in a folder, to make them programmatically accessible the `Nextflow` workflow, we will create a tar and zipped folder.

### Install [pigz](http://zlib.net/pigz/)

We will use [pigz](http://zlib.net/pigz/) after using tar.   

Searching [Anaconda](https://anaconda.org/conda-forge/pigz), we find a `conda install`.   Committed to managing our environments with `conda`.   A faster way of installing packages is to use a tool called `mamba` which we choose now to install.

```bash
conda install -c conda-forge mamba -y
```
Now we install [pigz](http://zlib.net/pigz/)

```bash
mamba install -c conda-forge pigz -y
```

### mass_spec folder

The fractionated proteomic data files were arranged in a folder in the AWS S3 bucket, each of these files begin with `120426_`.  We create a subdirectory, `mass_spec`, and then we tar that directory and then use  [pigz](http://zlib.net/pigz/) to zip it up.

Assuming we are in `Long-Read-Proteogenomics/data` directory.

#### i. Make the directory within `Long-Read-Proteogenomics/data` directory

```bash
mkdir mass_spec
```

#### ii. Move the mass_spec files to the directory just created

```bash
mv 120426_* mass_spec
```

#### iii. Inspect the files

Check to make sure the files have been successfully transferred.

#### iv. `tar` the files 

Using the unix command `tar`, we create 'c', verify 'v' the file
following specified with the option 'f'.  Assembled into the single command `tar cvf`.

The input to the `tar cvf` function is the directory `mass_spec`.
The output of the function is the file `mass_spec.tar`.

```bash
tar cvf mass_spec.tar mass_spec
```

#### v. [pigz](http://zlib.net/pigz/) the file for upload

The [pigz](http://zlib.net/pigz/) is a compression function that is rapid because it runs the compression in parallel. 
The input to this function is the `mass_spec.tar` file. 
The output of this function is the smaller `mass_spec.tar.gz`

```bash
pigz mass_spec.tar
```

### star_genome folder

The genome information required for `star` are in a folder in the AWS S3 Bucket.  To make this programmatically accessible through a `Channel` in a `Nextflow` workflow, we will use `tar` and `pigz` to create a `tar.gz`.

Assuming we are in `Long-Read-Proteogenomics/data` directory.

#### i. Make a directory within `Long-Read-Proteogenomics/data` directory

```bash
mkdir star_genome
```

#### ii. Move the `star_genome` files into the directory just created

Move the star_genome files from the `Long-Read-Proteogenomics/data` directory to the 
newly created directory `Long-Read-Proteogenomics/data/star_genome/`.

```bash
bash ../star_genome.sh
```

#### iii. Inspect the directory regarding its contents.

```bash
ls star_genome
star_genome/chrLength.txt
star_genome/chrNameLength.txt
star_genome/chrStart.txt
star_genome/geneInfo.tab
star_genome/exonGeTrInfo.tab
star_genome/exonInfo.tab
star_genome/sjdbInfo.txt
star_genome/genomeParameters.txt
star_genome/Log.out
star_genome/SA
star_genome/SAindex
star_genome/sjdbList.out.tab
star_genome/sjdbList.fromGTF.out.tab
star_genome/transcriptInfo.tab
```

#### iv. tar these files in the subdirectory

Using the tar command we create 'c', verify 'v' the file
following specified with the option 'f'.  Assembled into the single command `tar cvf`.
The input to the `tar cvf` function is the directory `star_genome`.
The output of the function is the file `star_genome.tar`.

```bash
tar cvf star_genome.tar star_genome
```

#### v. [pigz](http://zlib.net/pigz/) the file for compression prior to upload

The [pigz](http://zlib.net/pigz/) is a compression function that is rapid because it runs the compression in parallel. 
The input to this `parallel` zip function `pzip` is the `star_genome.tar` file.
The output of the [pigz](http://zlib.net/pigz/) function is the `star_genome.tar.gz` file.

```bash
pigz star_genome.tar
```

#### vi.  Inspect the contents (optional)

To inspect the contents of this now tar'd and zip'd file, run the following command.  The large `SA` file will take a while for the inspection to work through, so do not be surprised by the delay.

```bash
tar --list star_genome.tar.gz
```

## upload to Zenodo

To upload files to Zenodo, we take advantage of the application [zenodo-upload](https://github.com/jhpoelen/zenodo-upload).  Following the directions from the site: creating a repository, obtaining a `Zenodo Token` and executing the script to upload the files to the repository.

#### i. Create a Repository on Zenodo.

Either create a new repository, or create a new version of your repository, to put your data into the repository.  Be sure you have created and saved but not published the repository on Zenodo.  This puts the data in a state ready to receive new data.  

Pressing `Save` places the repository into a ready-state able to receive your data.  
Pressing `Publish`, places the version of the repository into a closed state, no longer ready to receive updates to the data repository.  

If you have pressed `Save` you are ready to go,
If you have pressed `Save` and `Publish` you are not and need to make a new version of the repository to put the repository back into the `ready-to-receive` state.

Press `Save` but do not `publish`. Now the repository is ready for data.

#### ii. Obtain a personal `Zenodo Token`.

Uploaded after reserving the `DOI` from Zenodo and getting a personal `Zenodo Token`, following the instructions [zenodo-upload](https://github.com/jhpoelen/zenodo-upload).  In short, under profile at the time of this writing on the upper right side of the browser window, scroll down and select `applications`.  There you will be able to generate a personal token.   The Deposit needs to be in a state where it can receive the data by this programmatic upload, so the new repository or version of the repository must be in a `Save` state and not yet `Published`.   The steps are outlined in further detail below.

Once you have your token, copy it, as it will disappear, the token is a secret key to allow you to upload the file.   Once copied, inside a terminal window, set an environment variable `ZENODO_TOKEN`.

```bash
export ZENODO_TOKEN=[`set to your own personal zenodo token`]
```

#### iii.  Upload the files to `Zenodo`

To run this pipeline using the `ZENODO` location does not require any downloading or movement of data.

Nevertheless, a bash script is provided to illustrate how one could pull the version of files that are stored within Zenodo.  This may be run from any unix terminal within an appropriately sized machine (2 vCPUs 16GB).

In general, it is not best practice to move data, download data or copy data.  Data movement should be minimized and discouraged to prevent integrity, limit costs, and ensoure consistency with the results.

The script documenting the uploading of the prepared files in `Long-Read-Proteogenomics/data` directory is `upload_to_zenodo.sh`
Run from `Long-Read-Proteogenomics` directory.

The pattern of upload is as follows:
```bash
cd zenodo-upload
./zenodo_upload.sh https://zenodo.org/deposit/5076056 ../Long-Read-Proteogenomics/data/gencode.v35.annotation.canonical.gtf
```

Since we have several files we are uploading, these are held within a simple bash script.

```bash
bash upload_to_zenodo.sh
 ```

## Download and Reconstruct

Now that the files are in Zenodo, we can download them and reconstruct the files as necessary

#### i.  Using Lifebit's CloudOS system, start a `jupyterlab notebook` (`2 CPUs with 16 GB RAM`), with 500 GB.

#### ii. Initialize `bash`

```bash
conda init bash
exec -l bash
```

#### iii. Create a `conda` environment

```bash
conda create -n lrp -y
conda activate lrp
```

#### iv. clone the `Long-Read-Proteogenomics` repository

```bash
git clone https://github.com/sheynkman-lab/Long-Read-Proteogenomics.git
cd Long-Read-Proteogenomics/data
```

#### v. Run the `download_and_reconstruct.sh` script (optional)

To run this pipeline using the `ZENODO` location does not require any downloading or movement of data.

Please run this script from the `Long-Read-Proteogenomics/data` directory.

```bash
bash download_and_reconstruct.sh
```
