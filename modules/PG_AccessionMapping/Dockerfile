# Full contents of Dockerfile

FROM continuumio/miniconda3:4.8.2
LABEL description="Base docker image with conda and util libraries"
ARG ENV_NAME="lrpg"

# Install mamba for faster installation in the subsequent step
# Install r-base for being able to run the install.R script
RUN conda install -c conda-forge mamba r-base -y

# Install the conda environment
COPY environment.yml /
RUN mamba env create --quiet --name ${ENV_NAME} --file /environment.yml && conda clean -a

# Install R packages that are possibly not available via conda
# COPY bin/install.R /
# RUN Rscript /install.R

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml

# Copy additional scripts from bin and add to PATH
RUN mkdir /opt/bin
COPY src/* /opt/bin/
RUN chmod +x /opt/bin/*
ENV PATH="$PATH:/opt/bin/"
