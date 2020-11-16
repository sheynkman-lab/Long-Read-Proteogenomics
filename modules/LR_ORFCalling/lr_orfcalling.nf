#!/usr/bin/env nextflow
/*
 * Copyright (c) 2020, Sheynkman Lab and the authors.
 *
 *   This file is part of 'proteogenomics-nf' a pipeline repository to run
 *   Gloria Sheynkman, Mike Shortreed and author's proteogenomics pipeline.
 *
 *   In this portion of the pipeline, LR_ORFCalling, the open reading frames are
 *   called, analyzed and filtered.
 *
 * @authors
 * Gloria Sheynkman
 * Ben Jordan
 * Anne Deslattes Mays (adeslat@scitechcon.org)
 */
# test
def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run lr_orfcalling.nf --fasta ../../data/jurkat_corrected.fasta -profile lr_orfcalling_nextflow.config
    
    Input files:
      --fasta                       path to the fasta file
      --run_name                    name for the run
 
    Transdecoder:                   no additional arguments required

    CPAT:                           
      --hexamer                     species specific hexamer file (cpat can also generate a new one)
      --gtf

    Other:
                                    (default: false)
      --max_cpus                    Maximum number of CPUs (int)
                                    (default: ?)  
      --max_memory                  Maximum memory (memory unit)
                                    (default: 80)
      --max_time                    Maximum time (time unit)
                                    (default: ?)
      --skipMultiQC                 Skip MultiQC (bool)
                                    (default: false)
      --outdir                      The output directory where the results will be saved (string)
                                    (default: directory where you submit the job)
      --mega_time                   Sets time limit for processes withLabel 'mega_memory' in the main.nf using the base.config (time unit)     
                                    (default: 20.h)
      --gc_disk_size                Only specific to google-cloud executor. Adds disk-space for few aggregative processes.
                                    (default: "200 GB" based on 100 samples. Simply add 2 x Number of Samples)

    See here for more info: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// main.nf

// get run name and date prefix for counts matrix and multiqc
run_name = params.run_name ? params.run_name + "_" : ""
date = new Date().format("MM-dd-yy")
run_prefix = run_name + date
outdir = run_name

log.info "lr_orfcalling - N F  ~  version 0.1"
log.info "====================================="
log.info "Run name                    : ${params.run_name}"
log.info "Date                        : ${date}"
log.info "Final prefix                : ${run_prefix}"
log.info "Fasta                       : ${params.fasta}"

if (params.fasta) {
  /*--------------------------------------------------
    TransDecoder for calling ORF on fasta file
  ---------------------------------------------------*/
  println "My fasta file is: ${params.fasta}"
  Channel
     .value(file(params.fasta))
     .ifEmpty { error "Cannot find any fasta file for parameter --fasta: ${params.fasta}" }
     .set { fasta }    

  process runTransDecoder {
    tag "runTransDecoder"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(fasta)

    output:
    file ("*.bed") into bed_channel
    file ("*.cds") into cds_channel
    file ("*.gff3") into gff3_channel
    file ("*.pep") into pep_channel

    script:
    """
    TransDecoder.LongOrfs -t ${params.fasta}
    TransDecoder.Predict -t ${params.fasta}    
    """
  }

}
