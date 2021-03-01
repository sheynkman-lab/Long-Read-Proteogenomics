## The folder `data`

### Why is this folder empty?

We use the folder **`data`** to store the artifacts that are generated during the analysis, so this folder should appear empty on GitHub. Large files cannot be pushed to GitHub as done with code but we keep this folder so that upon i) `git clone` of this repo and ii) acquiring the input files provided with the publication you will be able to reuse the Jupyter Notebooks without the need to change filepaths in the code.

### How to download data from Zenodo

Run these commands to automatically download data files:
`pip3 install zenodo_get`
`zenodo_get -w zenodo_data_links.txt <input here the latest version, see below*>`
`wget -i zenodo_data_links.txt`