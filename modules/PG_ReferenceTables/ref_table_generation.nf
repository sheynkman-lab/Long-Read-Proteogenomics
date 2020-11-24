#!/usr/bin/env nextflow
/*
 * Copyright (c) 2020, Sheynkman Lab and the authors.
 *
 *   This file is part of 'proteogenomics-nf' a pipeline repository to run
 *   Gloria Sheynkman, Mike Shortreed and author's proteogenomics pipeline.
 *
 *   In this portion of the pipeline, PG_ReferenceTables, helper tables of gene/isoform mapping and their length stats are generated. 
 *
 * @authors
 * Simi Kaur
 * Gloria Sheynkman
 */

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    // TODO - fill
    // example: nextflow run refined_db_generation.nf --orfs orf-testset-fraciton16.csv --seq jurkat_corrected.fasta --sample 'jurkat'
    
    Input files:
      --orfs          path to ORF coordinate input file location
      --seq           path to PacBio fasta sequence input file location
      --redundant     name of output redundant accession file location   
      --ctsv          name of output combined tsv file location
      --cfasta        ...
      --aggtsv 
      --aggfasta 
      --sample        name of unique identifier for the biosample
      
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
log.info "orfs : ${params.orfs}"
log.info "seq  : ${params.seq}"



  /*--------------------------------------------------
    Refined Protein Database Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.orfs))
     .ifEmpty { error "Cannot find orfs file for parameter --orfs: ${params.orfs}" }
     .set { ch_orfs }   
     
  Channel
     .value(file(params.seq))
     .ifEmpty { error "Cannot find any seq file for parameter --seq: ${params.seq}" }
     .set { ch_seq }  

  process generate_refined_database {
    tag "${orfs}, ${seq}"

    publishDir "${params.outdir}/refined_database/", mode: 'copy'

    input:
    file(orfs) from ch_orfs
    file(seq) from ch_seq
    
    output:
    file("*")
    
    script:
    """
    refine_orf.py \
    --orfs $orfs \
    --seq $seq \
    --redundant ${params.sample}_redundant_accessions.txt \
    --ctsv ${params.sample}_orf_combined.tsv \
    --cfasta ${params.sample}_orf_combined.fasta \
    --aggtsv ${params.sample}_orf_aggregated.tsv \
    --aggfasta ${params.sample}_orf_aggregated.fasta
    """
  }
