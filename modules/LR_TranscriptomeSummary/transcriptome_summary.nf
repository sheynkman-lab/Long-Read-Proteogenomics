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

log.info "transcriptome summary - N F  ~  version 0.1"
log.info "====================================="


sqanti_classification = false
  tpm = false
  ribo = false
  ensg_to_gene = false
  enst_to_isoname = false
  len_stats = false


  /*--------------------------------------------------
    Reference Table Generation 
  ---------------------------------------------------*/
  Channel
     .value(file(params.sqanti_classification))
     .ifEmpty { error "Cannot find gtf file for parameter --sqanti_classification: ${params.sqanti_classification}" }
     .set { ch_sqanti_classification }   
     
  Channel
     .value(file(params.tpm))
     .ifEmpty { error "Cannot find any tpm file for parameter --tpm: ${params.tpm}" }
     .set { ch_tpm }  
 
  Channel
     .value(file(params.ribo))
     .ifEmpty { error "Cannot find any file for parameter --ribo: ${params.ribo}" }
     .set { ch_ribo }  
     
  Channel
     .value(file(params.ensg_to_gene))
     .ifEmpty { error "Cannot find any file for parameter --ensg_to_gene: ${params.ensg_to_gene}" }
     .set { ch_ensg_to_gene } 
 
   Channel
     .value(file(params.enst_to_isoname))
     .ifEmpty { error "Cannot find any file for parameter --enst_to_isoname: ${params.enst_to_isoname}" }
     .set { ch_enst_to_isoname }  

   Channel
     .value(file(params.len_stats))
     .ifEmpty { error "Cannot find any file for parameter --len_stats: ${params.len_stats}" }
     .set { ch_len_stats }  
     
     

  process transcriptome_summary {

    publishDir "${params.outdir}/transcriptome_summary/", mode: 'copy'

    input:
    file(sqanti_classification) from ch_sqanti_classification
    file(tpm) from ch_tpm
    file(ribo) from ch_ribo
    file(ensg_to_gene) from ch_ensg_to_gene
    file(enst_to_isoname) from ch_enst_to_isoname
    file(len_stats) from ch_len_stats
    
    
    output:
    file("gene_level_tab.tsv") into ch_gene_level
    file("sqanti_isoform_info.tsv") into ch_sqanti_isoform_info
    file("*")
    
    script:
    """
    transcriptome_summary.py --sq_out $sqanti_classification --tpm $tpm --ribo $ribo --ensg_to_gene $ensg_to_gene --enst_to_isoname $enst_to_isoname --len_stats $len_stats 
    """
  }
