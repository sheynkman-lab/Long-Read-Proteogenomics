FROM continuumio/miniconda3:4.8.2

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y

#install transdecoder
RUN conda install -c bioconda transdecoder
