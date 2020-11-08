FROM continuumio/miniconda

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

#install cpat
RUN conda install -c bioconda cpat
