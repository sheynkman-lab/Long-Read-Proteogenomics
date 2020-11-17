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
    nextflow run lr_orfcalling.nf --fasta ../../data/jurkat_corrected.fasta -profile lr_orfcalling_nextflow.config
    
    Input files:
      --fasta                       path to the fasta file
      --trans_decoder               Boolean, defaults to true to execute the TransDecoder process
      
    Transdecoder:                   no additional arguments required

    Other:

      --max_cpus                    Maximum number of CPUs (int)
      --max_memory                  Maximum memory (memory unit)
      --max_time                    Maximum time (time unit)

    See here for more info: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

log.info "lr_orfcalling - N F  ~  version 0.1"
log.info "====================================="
log.info "Fasta                       : ${params.fasta}"



if (params.fasta && params.trans_decoder) {
  /*--------------------------------------------------
    TransDecoder for calling ORF on fasta file
  ---------------------------------------------------*/
  println "My fasta file is: ${params.fasta}"
  Channel
     .value(file(params.fasta))
     .ifEmpty { error "Cannot find any fasta file for parameter --fasta: ${params.fasta}" }
     .set { ch_fasta }    

  process run_transdecoder {
    tag "${fasta}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(fasta) from ch_fasta

    output:
    // alternative syntax to save all files in one channel:
    set file("*.bed"), file("*.cds"),  file("*.gff3"), file ("*.pep") into ch_print_to_check

    script:
    """
    TransDecoder.LongOrfs -t $fasta
    TransDecoder.Predict -t $fasta    
    """
  }
}


  process ch_print_to_check {
    echo true

    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file(bed), file(cds),  file(gff3), file(pep) from ch_print_to_check

    script:
    """
    ls -l 
    """
  }