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

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fasta my.fasta --gtf genome.gtf -profile base
    
    Input files:
      --fasta                       path to the fasta file
 
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
params.reads = false


// get run name and date prefix for counts matrix and multiqc
run_name = params.run_name ? params.run_name + "_" : ""
date = new Date().format("MM-dd-yy")
run_prefix = run_name + date

log.info "lr_orfcalling - N F  ~  version 0.1"
log.info "====================================="
log.info "Run name                    : ${params.run_name}"
log.info "Date                        : ${date}"
log.info "Final prefix                : ${run_prefix}"
log.info "Fasta                       : ${params.fasta}"
log.info "Outdir                      : ${params.outdir}"
log.info "Max CPUs                    : ${params.max_cpus}"
log.info "Max memory                  : ${params.max_memory}"
log.info "Max time                    : ${params.max_time}"
log.info "Mega time                   : ${params.mega_time}"
log.info "Google Cloud disk-space     : ${params.gc_disk_size}"

if (params.fasta) {
  /*--------------------------------------------------
    TransDecoder for calling ORF on fasta file
  ---------------------------------------------------*/
  println "My reads: ${params.reads}"

  process runTransDecoder {
    tag "$name"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(fasta) from ${params.fasta}

    output:
    file ${params.fasta}".{longest_orfs.pep,longest.gff3,longest_orfs.cds,longest_orfs.cds.top_500_longest}" into transdecoder_results 

    script:
    """
    TransDecoder.LongOrfs -t ${params.fasta}
    TransDecoder.Predict -t ${params.fasta}    
    """
  }

}