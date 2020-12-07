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
    // example: nextflow run ref_table_generation.nf --gtf gencode.gtf --fa gencode.fasta
    
    Input files:
      --gtf           path to the gencode gtf file
      --fa           path to the gencode transcript fasta file

    See here for more info: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

log.info "reference_tables - N F  ~  version 0.1"
log.info "====================================="
log.info "gtf : ${params.gtf}"
log.info "fa  : ${params.fa}"



  /*--------------------------------------------------
    Reference Table Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.gtf))
     .ifEmpty { error "Cannot find gtf file for parameter --gtf: ${params.gtf}" }
     .set { ch_gtf }   
     
  Channel
     .value(file(params.fa))
     .ifEmpty { error "Cannot find any fasta file for parameter --fa: ${params.fa}" }
     .set { ch_fa }  

  process generate_reference_tables {
    tag "${gtf}, ${fa}"

    publishDir "${params.outdir}/PG_ReferenceTables/", mode: 'copy'

    input:
    file(gtf) from ch_gtf
    file(fa) from ch_fa
    
    output:
    file("ensg_gene.tsv") into ch_ensg_gene
    file("enst_isoname.tsv") into ch_enst_isoname
    file("gene_ensp.tsv") into ch_gene_ensp
    file("gene_isoname.tsv") into ch_gene_isoname
    file("isoname_lens.tsv") into ch_isoname_lens
    file("gene_lens.tsv") into ch_gene_lens
    
    script:
    """
    python3 /opt/bin/prepare_reference_table.py \
        --gtf $gtf \
        --fa $fa \
        --ensg_gene ensg_gene.tsv \
        --enst_isoname enst_isoname.tsv \
        --gene_ensp gene_ensp.tsv \
        --gene_isoname gene_isoname.tsv \
        --isoname_lens isoname_lens.tsv \
        --gen_lens gene_lens.tsv
    """
  }
