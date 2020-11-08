# sheynkman-lab/Long-Read-Proteogenomics: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving sheynkman-lab/Long-Read-Proteogenomics.

We try to manage the required tasks for sheynkman-lab/Long-Read-Proteogenomics using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

> If you need help using or modifying sheynkman-lab/Long-Read-Proteogenomics then the best way to do so is by opening an issue in the repo with the label `help-needed`.

## Contribution workflow

If you'd like to write some code for sheynkman-lab/Long-Read-Proteogenomics, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [sheynkman-lab/Long-Read-Proteogenomics issues](https://github.com/sheynkman-lab/Long-Read-Proteogenomics/issues) to avoid duplicating work
    * If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [sheynkman-lab/Long-Read-Proteogenomics repository](https://github.com/sheynkman-lab/Long-Read-Proteogenomics) to your GitHub account
3. Make the necessary changes / additions within your forked repository
4. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

### Pipeline Tests

Each `nf-core` pipeline should be set up with a minimal set of test-data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

## Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

* On your own fork, make a new branch `patch` based on `upstream/master`.
* Fix the bug, and bump version (X.Y.Z+1).
* A PR should be made on `master` from patch to directly this particular bug.