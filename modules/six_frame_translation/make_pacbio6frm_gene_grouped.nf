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
    nextflow run make_pacbio6frm_gene_grouped.nf \
    --iso_annot SQANTI3_out/jurkat_classification.txt \
    --ensg_gene ensg_gene.tsv \
    --sample_fasta jurkat_corrected.fasta \
    --name jurkat
    --outdir ../


    Input files:
    --iso_annot       SQANTI3 jurkat classification file location
    --ensg_gene       ENSG_GENE file location from ReferenceTables
    --sample_fasta    Sample corrected fasta file

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

log.info "6frame - N F  ~  version 0.1"
log.info "====================================="
log.info "iso_annot : ${params.iso_annot}"
log.info "ensg_gene  : ${params.ensg_gene}"
log.info "sample_fasta : ${params.sample_fasta}"



  /*--------------------------------------------------
    Make PacBio 6-frame gene grouped
  ---------------------------------------------------*/
  Channel
     .value(file(params.iso_annot))
     .ifEmpty { error "Cannot find orfs file for parameter --iso_annot: ${params.iso_annot}" }
     .set { ch_iso_annot }   
     
  Channel
     .value(file(params.ensg_gene))
     .ifEmpty { error "Cannot find any seq file for parameter --ensg_gene: ${params.ensg_gene}" }
     .set { ch_ensg_gene } 

  Channel
     .value(file(params.sample_fasta))
     .ifEmpty { error "Cannot find any seq file for parameter --sample_fasta: ${params.sample_fasta}" }
     .set { ch_sample_fasta }  

  
  process make_pacbio_6frm_gene_grouped {
    tag "${iso_annot}, ${ensg_gene}"
    publishDir "${params.outdir}/pacbio_6frm_gene_grouped/", mode: 'copy'

    input:
    file(iso_annot) from ch_iso_annot
    file(ensg_gene) from ch_ensg_gene
    file(sample_fasta) from ch_sample_fasta
    
    output:
    file("${params.name}.6frame.fasta") into ch_6frm

    script:
    """
    make_pacbio6frm_gene_grouped.py --iso_annot $iso_annot --ensg_gene $ensg_gene --sample_fasta $sample_fasta --output_fasta ${params.name}.6frame.fasta
    """

  }

