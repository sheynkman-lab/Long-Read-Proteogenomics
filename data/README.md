## The folder `data`

### Why is this folder empty?

We use the folder **`data`** to store the artifacts that are generated during the analysis, so this folder should appear empty on GitHub. Large files cannot be pushed to GitHub as done with code but we keep this folder so that upon i) `git clone` of this repo and ii) acquiring the input files provided with the publication you will be able to reuse the Jupyter Notebooks without the need to change filepaths in the code.

### Results from the Manuscript
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5076056.svg)](https://doi.org/10.5281/zenodo.5076056)

For the publication, the results from the manuscript are found in the `9.1GB` file `LRPG-Manuscript-Results.tar.gz` in the Zenodo location.

### How to download data from Zenodo

As a practice, we encourage to always work within a controlled conda environment.

### Create a controlled environment

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
### Download the data from Zenodo

Within the `lrp` conda environment, install the `zenodo-get` package

```bash
conda install -c conda-forge zenodo_get
```

Test the installation:
```bash
(lrp) annedeslattesmays@Annes-MacBook-Pro data % zenodo_get -h
Usage: zenodo_get [options] RECORD_OR_DOI

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -c, --cite            print citation information
  -r RECORD, --record=RECORD
                        Zenodo record ID
  -d DOI, --doi=DOI     Zenodo DOI
  -m, --md5             Create md5sums.txt for verification.
  -w WGET, --wget=WGET  Create URL list for download managers. (Files will not
                        be downloaded.)
  -e, --continue-on-error
                        Continue with next file if error happens.
  -k, --keep            Keep files with invalid checksum. (Default: delete
                        them.)
  -n, --do-not-continue
                        Do not continue previous download attempt. (Default:
                        continue.)
  -R RETRY, --retry=RETRY
                        Retry on error N more times.
  -p PAUSE, --pause=PAUSE
                        Wait N second before retry attempt, e.g. 0.5
  -t TIMEOUT, --time-out=TIMEOUT
                        Set connection time-out. Default: 15 [sec].
  -o OUTDIR, --output-dir=OUTDIR
                        Output directory, created if necessary. Default:
                        current directory.
  -s, --sandbox         Use Zenodo Sandbox URL.
(lrp) annedeslattesmays@Annes-MacBook-Pro data % 
```

`zenodo_get` will be used to create the URL list for downloading or for download managers.
```bash
zenodo_get -w zenodo_data_links.txt
```

To get all of the links we use `wget`.
If you don't have `wget` you can use [Anaconda.org](https://anaconda.org) to get this utility
```bash
conda install -c conda-forge wget
```

Now we can get all the files that are on Zenodo
```bash
wget -i zenodo_data_links.txt
```


