# nf-core style template for Dockerfile
FROM continuumio/miniconda3:4.8.2
LABEL description="Base docker image containing util software requirements for the sheynkman-lab/Long-Read-Proteogenomics pipeline"

# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
RUN apt-get update \
 && apt-get install -y procps \
 && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/cpat/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name cpat > cpat.yml

