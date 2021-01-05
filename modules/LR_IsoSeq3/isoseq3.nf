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
        nextflow run isoseq3.nf \
        --css_reads ../css_reads.css \
        --gencode_fasta ../gencode_fasta.fasta \
        --primers_fasta ../primers_fasta.fasta 
    
    Input files:
    css_reads
    gencode_fasta
    primers_fasta
      

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
log.info "css_reads      : ${params.css_reads}"
log.info "gencode_fasta  : ${params.gencode_fasta}"
log.info "primers_fasta  : ${params.primers_fasta}"



  /*--------------------------------------------------
    Refined Protein Database Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.css_reads))
     .ifEmpty { error "Cannot find file for parameter --css_reads: ${params.css_reads}" }
     .set { ch_css_reads }   
     
  Channel
     .value(file(params.gencode_fasta))
     .ifEmpty { error "Cannot find any seq file for parameter --gencode_fasta: ${params.gencode_fasta}" }
     .set { ch_gencode_fasta }  

  Channel
     .value(file(params.primers_fasta))
     .ifEmpty { error "Cannot find any seq file for parameter --primers_fasta: ${params.primers_fasta}" }
     .set { ch_primers_fasta } 
  


  process isoseq3 {
    tag "${css_reads}, ${gencode_fasta}, ${primers_fasta}"

    publishDir "${params.outdir}/isoseq3/", mode: 'copy'

    input:
    file(css_reads) from ch_css_reads
    file(gencode_fasta) from ch_gencode_fasta
    file(primers_fasta) from ch_primers_fasta
    
    output:
    file("*")
    
    script:
    """
    # create an index
    pbindex $css_reads
 
    # module load isoseqenv
    lima --isoseq --dump-clips --peek-guess -j ${params.max_cpus} $css_reads $primers_fasta ${params.name}.demult.bam
    isoseq3 refine --require-polya ${params.name}.demult.NEB_5p--NEB_3p.subreadset.xml $primers_fasta ${params.name}.flnc.bam

    # clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
    isoseq3 cluster ${params.name}.flnc.bam ${params.name}.polished.bam --verbose --use-qvs

    # align reads to the genome, takes few minutes (40 core machine)
    pbmm2 align $gencode_fasta ${params.name}.polished.transcriptset.xml ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${params.max_cpus} --log-level INFO

    # collapse redundant reads
    isoseq3 collapse ${params.name}.aligned.bam ${params.name}.collapsed.gff
    """
  }

