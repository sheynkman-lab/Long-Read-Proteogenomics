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
    nextflow run refined_db_generation.nf --orfs orf-testset-fraciton16.csv --seq jurkat_corrected.fasta --sample 'jurkat'
    
    Input files:
      --gencode_fasta          path to gencode fasta file

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

log.info "pg_gencode_db - N F  ~  version 0.1"
log.info "====================================="
log.info "gencode_fasta : ${params.gencode_fasta}"



  /*--------------------------------------------------
    Refined Protein Database Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.gencode_fasta))
     .ifEmpty { error "Cannot find gencode fasta file for parameter --gencode_fasta: ${params.gencode_fasta}" }
     .set { ch_gencode_fasta }   
  
  
  process make_gencode_database {
    tag "${gencode_fasta}"

    publishDir "${params.outdir}/gencode_db/", mode: 'copy'

    input:
    file(gencode_fasta) from ch_gencode_fasta
    
    output:
    file("gencode.fasta") into ch_gencode_fasta_single
    file("gencode_isoname_clusters.tsv") into ch_gencode_isoname_clusters
    
    script:
    """
    make_gencode_database.py \
    --gencode_fasta $gencode_fasta \
    --output_fasta gencode.fasta \
    --output_cluster gencode_isoname_clusters.tsv
    """
  }

