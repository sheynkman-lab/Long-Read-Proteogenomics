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

    The typical command for running the orf module is as follows:
    nextflow run lr_orfcalling.nf --fasta https://zenodo.org/record/4278034/files/toy_for_christina.fasta.txt
    
    Input files:
      --fasta                [file] Path to the fasta file required for TransDecoder
      --trans_decoder        [bool] Boolean, defaults to true to execute the TransDecoder process


    See here for more info: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

log.info "--------------------------------------------------------------------------------------------"

def summary = [:]
if (params.fasta) summary['Fasta'] = "${params.fasta}"
summary['TransDecoder'] = "${params.trans_decoder}"


log.info summary.collect { k,v -> "${k.padRight( 18)}: $v" }.join("\n")
log.info "--------------------------------------------------------------------------------------------"

// Initialize channels based on parameters

// Stop early if the trans_decoder parameter is set to false
if (!params.trans_decoder) { exit 0, "Nothing to execute, set the --trans_decoder parameter to true if you wish to run this module."}

// Stop early if the fasta file (required for this process) is not provided
if (!params.fasta) { exit 1, "No fasta file found at the location ${params.fasta}. Please make sure the path to the file exists."}

// Create the channel for the fasta file if provided
if (params.fasta)  {ch_fasta = Channel.value(file(params.fasta, checkIfExists: true)) }

// Execute TransDecoder only when a --fasta file has been provided and --trans_decoder true
if (params.fasta && params.trans_decoder) {
  /*--------------------------------------------------
    TransDecoder for calling ORF on fasta file
  ---------------------------------------------------*/

  process run_transdecoder {
    tag "${fasta}"

    publishDir "${params.outdir}/transdecoder/", mode: 'copy'

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


  process print_to_check {
    echo true

    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file(bed), file(cds),  file(gff3), file(pep) from ch_print_to_check

    script:
    """
    echo "process 'print_to_check' staged input files, generated as outputs of process run_transdecoder:"
    ls -L
    """
  }

