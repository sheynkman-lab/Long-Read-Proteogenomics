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
nextflow run orf_calling.nf \
--orf_coord /mnt/shared/ubuntu/session_data/data/test_data/jurkat_cpat.ORF_prob_testset_fraction16.tsv \
--gencode_gtf /mnt/shared/ubuntu/session_data/data/test_data/gencode.v35.annotation.gtf \
--sample_gtf /mnt/shared/ubuntu/session_data/data/test_data/jurkat_corrected.gtf \
--pb_gene /mnt/shared/ubuntu/session_data/data/test_data/pb_to_gene.tsv \
--classification /mnt/shared/ubuntu/session_data/data/test_data/jurkat_classification.txt \
--sample_fasta /mnt/shared/ubuntu/session_data/data/test_data/jurkat_corrected.fasta \
--name jurkat_testset_fraction16
    
    Input files:
      | argument          | description                       | input module        |
      |-------------------|-----------------------------------|---------------------|
      | --orf_coord       | CPAT ORFs (TSV)                   | LR_CPAT             |
      | --gencode         | Gencode annotation (GTF)          |                     |
      | --sample_gtf      | SQANTI PB annotation (GTF)        | SQANTI              |
      | --pb_gene         | PacBio : Gencode Cross Ref       | PG_ReferenceTables  |
      | --classification  | SQANTI PB isoform classification  | SQANTI              |
      | --sample_fasta    | SQANTI PB sequence (FASTA)        | SQANTI              |
      

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
log.info "orf_coord      : ${params.orf_coord}"
log.info "gencode_gtf    : ${params.gencode_gtf}"
log.info "sample_gtf     : ${params.sample_gtf}"
log.info "pb_gene        : ${params.pb_gene}"
log.info "classification : ${params.classification}"
log.info "sample_fasta   : ${params.sample_fasta}"



  /*--------------------------------------------------
    Refined Protein Database Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.orf_coord))
     .ifEmpty { error "Cannot find orfs file for parameter --orf_coord: ${params.orf_coord}" }
     .set { ch_orf_coord }   
     
  Channel
     .value(file(params.gencode_gtf))
     .ifEmpty { error "Cannot find any seq file for parameter --gencode_gtf: ${params.gencode_gtf}" }
     .set { ch_gencode_gtf }  

  Channel
     .value(file(params.sample_gtf))
     .ifEmpty { error "Cannot find any seq file for parameter --sample_gtf: ${params.sample_gtf}" }
     .set { ch_sample_gtf } 
  
  Channel
     .value(file(params.pb_gene))
     .ifEmpty { error "Cannot find any seq file for parameter --seq: ${params.pb_gene}" }
     .set { ch_pb_gene } 
    
  Channel
     .value(file(params.classification))
     .ifEmpty { error "Cannot find any seq file for parameter --seq: ${params.classification}" }
     .set { ch_classification } 
  
  Channel
     .value(file(params.sample_fasta))
     .ifEmpty { error "Cannot find any seq file for parameter --seq: ${params.sample_fasta}" }
     .set { ch_sample_fasta } 
  


  process orf_calling {
    tag "${orf_coord}, ${gencode_gtf}, ${sample_gtf}, ${pb_gene}, ${classification}, ${sample_fasta} "

    publishDir "${params.outdir}/orf_calling/", mode: 'copy'

    input:
    file(orf_coord) from ch_orf_coord
    file(gencode_gtf) from ch_gencode_gtf
    file(sample_gtf) from ch_sample_gtf
    file(pb_gene) from ch_pb_gene
    file(classification) from ch_classification
    file(sample_fasta) from ch_sample_fasta
    
    output:
    file("*")
    
    script:
    """
    orf_calling.py --orf_coord $orf_coord --gencode $gencode_gtf --sample_gtf $sample_gtf --pb_gene $pb_gene --classification $classification --sample_fasta $sample_fasta --output ${params.name}_orf_called.tsv
    """
  }

