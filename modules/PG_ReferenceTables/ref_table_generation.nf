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
    // example: nextflow run ref_table_generation.nf --gencode_gtf gencode.gtf --gencode_fasta gencode.fasta
    
    Input files:
      --gencode_gtf           path to the gencode gtf file
      --gencode_fasta           path to the gencode transcript fasta file

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
log.info "gencode_gtf : ${params.gencode_gtf}"
log.info "gencode_fasta  : ${params.gencode_fasta}"



  /*--------------------------------------------------
    Reference Table Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.gencode_gtf))
     .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.gencode_gtf}" }
     .set { ch_gencode_gtf }   
     
  Channel
     .value(file(params.gencode_fasta))
     .ifEmpty { error "Cannot find any gencode_fasta file for parameter --gencode_fasta: ${params.gencode_fasta}" }
     .set { ch_gencode_fasta }  

  process generate_reference_tables {
    tag "${gencode_gtf}, ${gencode_fasta}"

    publishDir "${params.outdir}/PG_ReferenceTables/", mode: 'copy'

    input:
    file(gencode_gtf) from ch_gencode_gtf
    file(gencode_fasta) from ch_gencode_fasta
    
    output:
    file("ensg_gene.tsv") into ch_ensg_gene
    file("enst_isoname.tsv") into ch_enst_isoname
    file("gene_ensp.tsv") into ch_gene_ensp
    file("gene_isoname.tsv") into ch_gene_isoname
    file("isoname_lens.tsv") into ch_isoname_lens
    file("gene_lens.tsv") into ch_gene_lens
    file("protein_coding_genes.txt") into ch_protein_coding_genes
    
    script:
    """
    prepare_reference_tables.py \
    --gtf $gencode_gtf \
    --fa $gencode_fasta \
    --ensg_gene ensg_gene.tsv \
    --enst_isoname enst_isoname.tsv \
    --gene_ensp gene_ensp.tsv \
    --gene_isoname gene_isoname.tsv \
    --isoname_lens isoname_lens.tsv \
    --gen_lens gene_lens.tsv \
    --protein_coding_genes protein_coding_genes.txt
    """
  }
