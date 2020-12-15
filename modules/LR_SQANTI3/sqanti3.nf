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
log.info "fl_count       : ${params.fl_count}"
log.info "gencode_gtf    : ${params.gencode_gtf}"
log.info "gencode_fasta  : ${params.gencode_fasta}"
log.info "sample_gtf     : ${params.sample_gtf}"




  /*--------------------------------------------------
    Refined Protein Database Generation 
  ---------------------------------------------------*/
  
  
  Channel
     .value(file(params.fl_count))
     .ifEmpty { error "Cannot find orfs file for parameter --fl_count: ${params.fl_count}" }
     .set { ch_fl_count }   
     
  Channel
     .value(file(params.gencode_gtf))
     .ifEmpty { error "Cannot find any seq file for parameter --gencode_gtf: ${params.gencode_gtf}" }
     .set { ch_gencode_gtf }  

    Channel
     .value(file(params.gencode_fasta))
     .ifEmpty { error "Cannot find any seq file for parameter --gencode_fasta: ${params.gencode_fasta}" }
     .set { ch_gencode_fasta }  


  Channel
     .value(file(params.sample_gtf))
     .ifEmpty { error "Cannot find any seq file for parameter --sample_gtf: ${params.sample_gtf}" }
     .set { ch_sample_gtf } 

  Channel
     .value(file(params.sample_fasta))
     .ifEmpty { error "Cannot find any fasta file for parameter --sample_fasta: ${params.sample_fasta}" }
     .set { ch_sample_fasta } 
  

  


  process sqanti3 {
    tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"

    publishDir "${params.outdir}/sqanti3/", mode: 'copy'

    input:
    
    file(fl_count) from ch_fl_count
    file(gencode_gtf) from ch_gencode_gtf
    file(gencode_fasta) from ch_gencode_fasta
    file(sample_fasta) from ch_sample_fasta
    
    
    output:
    file("*")
    
    script:
    """
    sqanti3_qc $sample_fasta $gencode_gtf $gencode_fasta -o ${params.name} -d SQANTI3_out/ --fl_count $fl_count

    """
    //
  }