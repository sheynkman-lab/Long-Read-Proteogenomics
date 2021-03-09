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