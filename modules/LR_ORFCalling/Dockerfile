# nf-core style template for Dockerfile
FROM continuumio/miniconda3:4.8.2
LABEL description="Basic docker image containing util software requirements for the sheynkman-lab/Long-Read-Proteogenomics pipeline ORF Calling Module"
ARG ENV_NAME="orf-calling"
# Install procps so that Nextflow can poll CPU usage and
# deep clean the apt cache to reduce image/layer size
# ! Leave this as is for any module
RUN apt-get update \
 && apt-get install -y procps \
 && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# ! Leave this as is for any module
# Install the conda environment
COPY environment.yml /
RUN conda env create --name ${ENV_NAME} --quiet --file /environment.yml && conda clean -a

# ! Update the envs/<name> based on the environment.yml file
# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ${ENV_NAME} > ${ENV_NAME}.yml


# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

# Copy additional scripts from bin and add to PATH
RUN mkdir /opt/bin
COPY bin/* /opt/bin/
RUN chmod +x /opt/bin/*
ENV PATH="$PATH:/opt/bin/"
