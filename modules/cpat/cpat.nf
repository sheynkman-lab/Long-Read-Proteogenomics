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
      | --pb_gene         | PacBio : Gencode Cross Ref        | PG_ReferenceTables  |
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

log.info "cpat - N F  ~  version 0.1"
log.info "====================================="
log.info "hexamer       : ${params.hexamer}"
log.info "logit model   : ${params.logit_model}"
log.info "sample fasta  : ${params.sample_fasta}"




  /*--------------------------------------------------
    CPAT
  ---------------------------------------------------*/
  
  
  Channel
     .value(file(params.hexamer))
     .ifEmpty { error "Cannot find headmer file for parameter --hexamer: ${params.hexamer}" }
     .set { ch_hexamer }   
     
  Channel
     .value(file(params.logit_model))
     .ifEmpty { error "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}" }
     .set { ch_logit_model }  

  Channel
     .value(file(params.sample_fasta))
     .ifEmpty { error "Cannot find any fasta file for parameter --sample_fasta: ${params.sample_fasta}" }
     .set { ch_sample_fasta } 
  

  


  process cpat {
    tag "${hexamer}, ${logit_model}, ${sample_fasta}"

    publishDir "${params.outdir}/cpat/", mode: 'copy'

    input:
    
    file(hexamer) from ch_hexamer
    file(logit_model) from ch_logit_model
    file(sample_fasta) from ch_sample_fasta
    
    
    output:
    file("${params.name}.ORF_prob.tsv") into ch_cpat_orfs
    file("*")
    
    script:
    """
    cpat.py \
    -x $hexamer \
    -d $logit_model \
    -g $sample_fasta \
    --min-orf=${params.min_orf} \
    --top-orf=50 \
    -o ${params.name} \
    1> ${params.name}_cpat.output \
    2> ${params.name}_cpat.error
    """
  }

