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
     .value(file(params.protein_coding_genes))
     .ifEmpty { error "Cannot find orfs file for parameter --orfs: ${params.orfs}" }
     .set { ch_protein_coding_genes }   
     
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
    --pb_fasta $sample_fasta \
    --redundant ${params.sample}_redundant_accessions.txt \
    --combined_tsv ${params.sample}_orf_combined.tsv \
    --combined_fasta ${params.sample}_orf_combined.fasta \
    --agg_tsv ${params.sample}_orf_aggregated.tsv \
    --agg_fasta ${params.sample}_orf_aggregated.fasta \
    --protein_coding_only ${params.protein_coding_only} \
    --protein_coding_genes $ch_protein_coding_genes \
    --cutoff ${params.refine_cutoff} \
    
    """
  }

