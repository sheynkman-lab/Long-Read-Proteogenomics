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
log.info "sample_gtf : ${params.sample_gtf}"
log.info "agg orfs : ${params.agg_orfs}"
log.info "refined orfs ${params.refined_orfs}"
log.info "pb_gene ${params.pb_gene}"



  /*--------------------------------------------------
    Refined Protein Database Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.sample_gtf))
     .ifEmpty { error "Cannot find file for parameter --sample_gtf: ${params.sample_gtf}" }
     .set { ch_sample_gtf }   

  Channel
     .value(file(params.agg_orfs))
     .ifEmpty { error "Cannot find file for parameter --agg_orfs: ${params.agg_orfs}" }
     .set { ch_agg_orfs }  
  Channel
     .value(file(params.refined_orfs))
     .ifEmpty { error "Cannot find file for parameter --refined_orfs: ${params.refined_orfs}" }
     .set { ch_refined_orfs}  
  
  Channel
     .value(file(params.pb_gene))
     .ifEmpty { error "Cannot find file for parameter --pb_gene: ${params.pb_gene}" }
     .set { ch_pb_gene }  
  
  
  process make_pacbio_cds_gtf {
    tag "${gencode_fasta}"

    publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

    input:
    file(sample_gtf) from ch_sample_gtf
    file(agg_orfs) from ch_agg_orfs
    file(refined_orfs) from ch_refined_orfs
    file(ch_pb_gene) from ch_pb_gene
    
    output:
    file("*")
    
    script:
    """
    make_pacbio_cds_gtf.py --sample_gtf $sample_gtf --agg_orfs $agg_orfs --refined_orfs $refined_orfs --pb_gene $pb_gene --output_cds ${params.name}_cds.gtf
    """
  }

