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

log.info "Longread Proteogenomics - N F  ~  version 0.1"
log.info "====================================="
log.info "orf_coord      : ${params.orf_coord}"
log.info "gencode_gtf    : ${params.gencode_gtf}"
log.info "sample_gtf     : ${params.sample_gtf}"
log.info "pb_gene        : ${params.pb_gene}"
log.info "classification : ${params.classification}"
log.info "sample_fasta   : ${params.sample_fasta}"




  Channel
     .value(file(params.cpat_orfs))
     .ifEmpty { error "Cannot find orfs file for parameter --cpat_orfs: ${params.cpat_orfs}" }
     .set { ch_cpat_orfs }   
     
  Channel
     .value(file(params.gencode_gtf))
     .ifEmpty { error "Cannot find any seq file for parameter --gencode_gtf: ${params.gencode_gtf}" }
     .set { ch_gencode_gtf }  

  Channel
     .value(file(params.sample_gtf))
     .ifEmpty { error "Cannot find any file for parameter --sample_gtf: ${params.sample_gtf}" }
     .set { ch_sample_gtf } 
     
  Channel
     .value(file(params.sample_fasta))
     .ifEmpty { error "Cannot find any file for parameter --sample_fasta: ${params.sample_fasta}" }
     .set { ch_sample_fasta } 
  
  Channel
     .value(file(params.pb_gene))
     .ifEmpty { error "Cannot find any file for parameter --pb_gene: ${params.pb_gene}" }
     .set { ch_pb_gene } 
    
  Channel
     .value(file(params.classification))
     .ifEmpty { error "Cannot find any file for parameter --classification: ${params.classification}" }
     .set { ch_classification } 
     
  Channel
     .value(file(params.protein_coding_genes))
     .ifEmpty { error "Cannot find any file for parameter --protein_coding_genes: ${params.protein_coding_genes}" }
     .set { ch_protein_coding_genes } 
  
/*--------------------------------------------------
ORF Calling 
---------------------------------------------------*/

  process orf_calling {
    tag "${orf_coord}, ${gencode_gtf}, ${sample_gtf}, ${pb_gene}, ${classification}, ${sample_fasta} "

    publishDir "${params.outdir}/orf_calling/", mode: 'copy'

    input:
    file(cpat_orfs) from ch_cpat_orfs
    file(gencode_gtf) from ch_gencode_gtf
    file(sample_gtf) from ch_sample_gtf
    file(sample_fasta) from ch_sample_fasta
    file(pb_gene) from ch_pb_gene
    file(classification) from ch_classification
    
    
    output:
    file("${params.name}_best_orf.tsv") into ch_best_orf
    
    script:
    """
    orf_calling.py --orf_coord $cpat_orfs --gencode $gencode_gtf --sample_gtf $sample_gtf --pb_gene $pb_gene --classification $classification --sample_fasta $sample_fasta --output ${params.name}_best_orf.tsv
    """
  }
  
/*--------------------------------------------------
Refined Protein Database Generation 
---------------------------------------------------*/
   
  process generate_refined_database {
    tag "${orfs}, ${seq}"

    publishDir "${params.outdir}/refined_database/", mode: 'copy'

    input:
    file(best_orfs) from ch_best_orf
    file(sample_fasta) from ch_sample_fasta
    file(protein_coding_genes) from ch_protein_coding_genes
    
    output:
    file("*")
    
    script:
    """
    refine_orf.py \
    --orfs $best_orfs \
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
