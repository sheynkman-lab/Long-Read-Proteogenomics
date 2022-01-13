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
      nextflow $NEXTFLOW_DIR/sqanti3.nf \
      --name $NAME \
      --genome_fasta $GENOME_FASTA \
      --gencode_gtf   $GENCODE_GTF \
      --sample_gtf    $SAMPLE_GTF \
      --fl_count      $FL_COUNT
    
    Input files:
      | argument          | description                       |
      |-------------------|-----------------------------------|
      | --fl_count       | isoseq full length counts               |
      | --gencode_gtf         | Gencode annotation (GTF)          |
      | --sample_gtf      | IsoSeq annotation (GTF)        |
      | --gencode_fasta         | GGENCODE   |
      | --sample_gtf  | SQANTI PB isoform classification  |
      | --sample_fasta    | SQANTI PB sequence (FASTA)        |
      

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
log.info "genome_fasta  : ${params.genome_fasta}"
log.info "sample_gtf     : ${params.sample_gtf}"



Channel
    .value(file(params.fl_count))
    .ifEmpty { error "Cannot find orfs file for parameter --fl_count: ${params.fl_count}" }
    .set { ch_fl_count }   
    
Channel
    .value(file(params.gencode_gtf))
    .ifEmpty { error "Cannot find any seq file for parameter --gencode_gtf: ${params.gencode_gtf}" }
    .set { ch_gencode_gtf }  

Channel
    .value(file(params.sample_gtf))
    .ifEmpty { error "Cannot find any seq file for parameter --sample_gtf: ${params.sample_gtf}" }
    .set { ch_sample_gtf } 


if (!params.genome_fasta) exit 1, "Cannot find any file for parameter --genome_fasta: ${params.genome_fasta}"
if (params.genome_fasta.endsWith('.gz')){
   ch_genome_fasta = Channel.value(file(params.genome_fasta))
} else {
   ch_genome_fasta_uncompressed = Channel.value(file(params.genome_fasta))
}



if (params.genome_fasta.endsWith('.gz')) {
   process gunzip_gencome_fasta {
   tag "decompress gzipped genome fasta"
   cpus 1

   input:
   file(genome_fasta) from ch_genome_fasta

   output:
   file("*.{fa,fasta}") into ch_genome_fasta_uncompressed

   script:
   """
   gunzip -f ${genome_fasta}
   """
   }
}



process sqanti3 {
  tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"

  publishDir "${params.outdir}/sqanti3/", mode: 'copy'

  input:
  
  file(fl_count) from ch_fl_count
  file(gencode_gtf) from ch_gencode_gtf
  file(genome_fasta) from ch_genome_fasta_uncompressed
  file(sample_fasta) from ch_sample_fasta
  
  
  output:
  file("*")
  
  script:
  """
    sqanti3_qc.py \
    $sample_gtf \
    $gencode_gtf \
    $genome_fasta \
    --skipORF \
    -o ${params.name} \
    --fl_count $fl_count  \
    --gtf
    """
}