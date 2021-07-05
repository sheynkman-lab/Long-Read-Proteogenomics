# Using Zenodo

To make the data more accessible and FAIR, the indexed files were transfered to Zenodo using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload) from the `University of Virginia's Gloria Sheynkman Lab` Amazon `S3` buckets.

The data were prepared and stored in the development of the `Sheynkman Lab Long-read Proteogenomics Pipeline`

## Prepare your environment

<p align="center"><img src="https://github.com/lifebit-ai/dry-bench-skills-for-researchers/blob/adds-mini-courses/assets/lifebitCloudOS.png"  width="250" align="right" ></p>

We used `Lifebit's CloudOS`, use conda to organize our environment, this maybe repeated in an appropriate sized terminal, running a unix environment.

i. initialize the bash environment

```bash
conda init bash
exec -l bash
```

ii. create and activate a new conda environment `lrp`.
```bash
conda create -n lrp
conda activate lrp
```

Now with the environment ready - begin the pre-requisites

## Clone the [zenodo-upload](https://github.com/jhpoelen/zenodo-upload) repository

```bash
git clone https://github.com/jhpoelen/zenodo-upload.git
```

## Satisfy the Pre-requisites for zenodo upload


i. Install [`jq`](https://stedolan.github.io/jq/).  

```bash
conda install -c conda-forge jq -y
```

ii. Installed `curl`

```bash
conda install -c conda-forge curl -y
```

iii. Bash

already satisfied.

## Preparing to use [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload)

## Transfer data from an S3 bucket

Data were on aws.  Add `awscli` the command line library for AWS S3 buckets.
For organizational principals, begin in the data directory of the recently cloned file.

i. install command-line interface to Amazon Web Services `awscli`

```bash
conda install -c conda-forge awscli -y
```

ii. clone the repository

```bash
git clone https://github.com/sheynkman-lab/Long-Read-Proteogenomics.git
```

iii. download the files

```bash
cd Long-Read-Proteogenomics/data
bash ../download_aws.sh
```

The download script has flattened the hierarchy of files.  Ensuring the scripts run appropriately, we need to recreate the directory structure where many of these files were previously arranged into folders.  These folders will be replaced by tar'd and zipped files that are stored in Zenodo. 

We do this with tar (system provided) and [pigz](http://zlib.net/pigz/)

## tar and pigz

For each previously arranged file that was in a folder, we will create a tar and zipped folder. We will use [pigz](http://zlib.net/pigz/) after using tar.   

Searching [Anaconda](https://anaconda.org/conda-forge/pigz), we find a `conda install`.   Committed to managing our environments with `conda`.   A faster way of installing packages is to use a tool called `mamba` which we choose now to install.

```bash
conda install -c conda-forge mamba -y
```
Now we install [pigz](http://zlib.net/pigz/)

```bash
mamba install -c conda-forge pigz -y
```

## mass_spec folder

The fractionated proteomic data files were arranged in a folder in the AWS S3 bucket, each of these files begin with `120426_`.  We create a subdirectory, `mass_spec`, and then we tar that directory and then use  [pigz](http://zlib.net/pigz/) to zip it up.

Assuming we are in `Long-Read-Proteogenomics/data` directory.

i. Make the directory within `Long-Read-Proteogenomics/data` directory

```bash
mkdir mass_spec
```

ii. Move the mass_spec files to the directory just created

```bash
mv 120426_* mass_spec
```

iii. Tar the files 

```bash
tar cvf mass_spec.tar mass_spec
```

iv. Pigz zip the tar file

```bash
pigz mass_spec.tar
```

## star_genome folder

The genome information required for `star` are in a folder in the AWS S3 Bucket.
Assuming we are in `Long-Read-Proteogenomics/data` directory.

i. Make a directory within `Long-Read-Proteogenomics/data` directory

```bash
mkdir star_genome
```

ii. Move the `star_genome` files into the directory just created

Move the star_genome files from the `Long-Read-Proteogenomics/data` directory to the 
newly created directory `Long-Read-Proteogenomics/data/star_genome/`.

```bash
bash ../star_genome.sh
```

iii. Inspect and then tar

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

iv. tar these files in the subdirectory

```bash
tar cvf star_genome.tar star_genome
```

v. pigz zip the file for upload

```bash
pigz star_genome.tar
```


## upload to Zenodo

Uploaded after reserving the `DOI` from Zenodo and getting a personal `zenodo token`, following the instructions [zenodo-upload](https://github.com/jhpoelen/zenodo-upload), I set the ZENODO_TOKEN environment variable.  

```bash
export ZENODO_TOKEN=[`set to your own personal zenodo token`]
```

## Create or make a new version of a Zenodo repository and obtain a personal upload token

First, be sure you have created and saved but not published the repository on Zenodo.  
If you need to make an update to the respository, you need to make a new version and repeat to `Save` but do not `publish`.   This allows you to upload the files.  Second, you will need to obtain a `Zenodo token`

Therefore two requirements need to be satisfied to upload data to Zenodo

i. You have an account on Zenodo and you have obtained a token.

Follow the steps for obtaining a token from [jhpoelen](https://github.com/jhpoelen/zenodo-upload#usage)

ii. Create the repository. upload to **Zenodo**

There are several files in this so we have a script, `upload-to-zenodo.sh`.   The `zenodo-upload` procedure has execution from the cloned repository of `zenodo-upload

A script was created to upload the files that have been downloaded and the two files that have been tar'd and zipped (`mass_spec.tar.gz` and `star_genome.tar.gz`)

Now we upload them to the new deposit, save and publish on the Zenodo site

i. From the `Long-Read-Proteogenomics` subdirectory

The script was designed to run from the `Long-Read-Proteogenomics` subdirectory, move there if not already there

```bash
bash upload_to_zenodo.sh
 ```

ii. 
## downloading from Zenodo

i.  Using Lifebit's CloudOS system, start a `jupyterlab notebook` (), with 1500 GB.

ii. Start a `bash command shell`

iii. Initialize `bash`

```bash
conda init bash
exec -l bash
```

iv. Create a `conda` environment

```bash
conda create -n lrp -y
conda activate lrp
```

v. clone the `Long-Read-Proteogenomics` repository

```bash
git clone https://github.com/sheynkman-lab/Long-Read-Proteogenomics.git
cd Long-Read-Proteogenomics/data
```

v. Download and Reconstruct

A bash script was created to pull the version of files that are stored within Zenodo


```bash
bash download_and_reconstruct.sh
```
