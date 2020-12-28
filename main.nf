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
// log.info "orf_coord      : ${params.orf_coord}"
// log.info "gencode_gtf    : ${params.gencode_gtf}"
// log.info "sample_gtf     : ${params.sample_gtf}"
// log.info "pb_gene        : ${params.pb_gene}"
// log.info "classification : ${params.classification}"
// log.info "sample_fasta   : ${params.sample_fasta}"




//   Channel
//      .value(file(params.cpat_orfs))
//      .ifEmpty { error "Cannot find orfs file for parameter --cpat_orfs: ${params.cpat_orfs}" }
//      .set { ch_cpat_orfs }   
     
//   Channel
//      .value(file(params.gencode_gtf))
//      .ifEmpty { error "Cannot find any seq file for parameter --gencode_gtf: ${params.gencode_gtf}" }
//      .set { ch_gencode_gtf }  

//   Channel
//      .value(file(params.sample_gtf))
//      .ifEmpty { error "Cannot find any file for parameter --sample_gtf: ${params.sample_gtf}" }
//      .set { ch_sample_gtf } 
     
//   Channel
//      .value(file(params.sample_fasta))
//      .ifEmpty { error "Cannot find any file for parameter --sample_fasta: ${params.sample_fasta}" }
//      .set { ch_sample_fasta } 
  
//   Channel
//      .value(file(params.pb_gene))
//      .ifEmpty { error "Cannot find any file for parameter --pb_gene: ${params.pb_gene}" }
//      .set { ch_pb_gene } 
    
//   Channel
//      .value(file(params.classification))
//      .ifEmpty { error "Cannot find any file for parameter --classification: ${params.classification}" }
//      .set { ch_classification } 
     
//   Channel
//      .value(file(params.protein_coding_genes))
//      .ifEmpty { error "Cannot find any file for parameter --protein_coding_genes: ${params.protein_coding_genes}" }
//      .set { ch_protein_coding_genes } 

Channel
    .value(file(params.gencode_gtf))
    .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.gencode_gtf}" }
    .set { ch_gencode_gtf }   
    
Channel
    .value(file(params.gencode_transcript_fasta))
    .ifEmpty { error "Cannot find any gencode_fasta file for parameter --gencode_fasta: ${params.gencode_transcript_fasta}" }
    .set { ch_gencode_transcript_fasta }  


Channel
    .value(file(params.sample_css))
    .ifEmpty { error "Cannot find file for parameter --sample_css: ${params.sample_css}" }
    .set { ch_sample_css }   
    
Channel
    .value(file(params.gencode_fasta))
    .ifEmpty { error "Cannot find any seq file for parameter --gencode_fasta: ${params.gencode_fasta}" }
    .set { ch_gencode_fasta }  

Channel
    .value(file(params.primers_fasta))
    .ifEmpty { error "Cannot find any seq file for parameter --primers_fasta: ${params.primers_fasta}" }
    .set { ch_primers_fasta } 

Channel
     .value(file(params.hexamer))
     .ifEmpty { error "Cannot find headmer file for parameter --hexamer: ${params.hexamer}" }
     .set { ch_hexamer }   
     
  Channel
     .value(file(params.logit_model))
     .ifEmpty { error "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}" }
     .set { ch_logit_model } 
/*--------------------------------------------------
Reference Tables 
---------------------------------------------------*/


process generate_reference_tables {
  tag "${gencode_gtf}, ${gencode_transcript_fasta}"

  publishDir "${params.outdir}/PG_ReferenceTables/", mode: 'copy'

  input:
  file(gencode_gtf) from ch_gencode_gtf
  file(gencode_transcript_fasta) from ch_gencode_transcript_fasta
  
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
  --fa $gencode_transcript_fasta \
  --ensg_gene ensg_gene.tsv \
  --enst_isoname enst_isoname.tsv \
  --gene_ensp gene_ensp.tsv \
  --gene_isoname gene_isoname.tsv \
  --isoname_lens isoname_lens.tsv \
  --gene_lens gene_lens.tsv \
  --protein_coding_genes protein_coding_genes.txt
  """
}

/*--------------------------------------------------
Gencode Database
---------------------------------------------------*/


/*--------------------------------------------------
Accession Mapping 
---------------------------------------------------*/

/*--------------------------------------------------
IsoSeq3
---------------------------------------------------*/
process isoseq3 {
  tag "${sample_css}, ${gencode_fasta}, ${primers_fasta}"

  publishDir "${params.outdir}/isoseq3/", mode: 'copy'

  input:
  file(sample_css) from ch_sample_css
  file(gencode_fasta) from ch_gencode_fasta
  file(primers_fasta) from ch_primers_fasta
  
  output:
  file("*")
  file("${params.name}.collapsed.gff") into ch_sample_gtf
  file("${params.name}.collapsed.fasta") into ch_sample_fasta
  file("${params.name}.collapsed.abundance.txt") into ch_fl_count
  script:
  """
  # create an index
  pbindex $sample_css

  # module load isoseqenv
  lima --isoseq --dump-clips --peek-guess -j ${params.max_cpus} $sample_css $primers_fasta ${params.name}.demult.bam
  isoseq3 refine --require-polya ${params.name}.demult.NEB_5p--NEB_3p.subreadset.xml $primers_fasta ${params.name}.flnc.bam

  # clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
  isoseq3 cluster ${params.name}.flnc.bam ${params.name}.polished.bam --verbose --use-qvs

  # align reads to the genome, takes few minutes (40 core machine)
  pbmm2 align $gencode_fasta ${params.name}.polished.transcriptset.xml ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${params.max_cpus} --log-level INFO

  # collapse redundant reads
  isoseq3 collapse ${params.name}.aligned.bam ${params.name}.collapsed.gff
  """
}
/*--------------------------------------------------
Six-Frame Translation
---------------------------------------------------*/

/*--------------------------------------------------
SQANTI3
---------------------------------------------------*/
  process sqanti3 {
    tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"

    publishDir "${params.outdir}/sqanti3/", mode: 'copy'

    input:
    
    file(fl_count) from ch_fl_count
    file(gencode_gtf) from ch_gencode_gtf
    file(gencode_fasta) from ch_gencode_fasta
    file(sample_gtf) from ch_sample_gtf
    
    
    output:
    file("*")
    
    script:
    """
    sqanti3_qc.py \
    $sample_gtf \
    $gencode_gtf \
    $gencode_fasta \
    --skipORF \
    -o ${params.name} \
    --fl_count $fl_count  \
    --gtf
    """
    //
  }
/*--------------------------------------------------
Transcriptome Summary 
---------------------------------------------------*/

/*--------------------------------------------------
CPAT
---------------------------------------------------*/
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
  --min-orf=50 \
  --top-orf=50 \
  -o ${params.name} \
  1> ${params.name}_cpat.output \
  2> ${params.name}_cpat.error

  """
}
/*--------------------------------------------------
ORF Calling 
---------------------------------------------------*/
// process orf_calling {
//   tag "${orf_coord}, ${gencode_gtf}, ${sample_gtf}, ${pb_gene}, ${classification}, ${sample_fasta} "

//   publishDir "${params.outdir}/orf_calling/", mode: 'copy'

//   input:
//   file(cpat_orfs) from ch_cpat_orfs
//   file(gencode_gtf) from ch_gencode_gtf
//   file(sample_gtf) from ch_sample_gtf
//   file(sample_fasta) from ch_sample_fasta
//   file(pb_gene) from ch_pb_gene
//   file(classification) from ch_classification
  
  
//   output:
//   file("${params.name}_best_orf.tsv") into ch_best_orf
  
//   script:
//   """
//   orf_calling.py \
//   --orf_coord $cpat_orfs \
//   --gencode $gencode_gtf \
//   --sample_gtf $sample_gtf \
//   --pb_gene $pb_gene \
//   --classification $classification \
//   --sample_fasta $sample_fasta \
//   --output ${params.name}_best_orf.tsv
//   """
// }

/*--------------------------------------------------
Refined DB Generation 
---------------------------------------------------*/
// process generate_refined_database {
//   tag "${orfs}, ${seq}"

//   publishDir "${params.outdir}/refined_database/", mode: 'copy'

//   input:
//   file(best_orfs) from ch_best_orf
//   file(sample_fasta) from ch_sample_fasta
//   file(protein_coding_genes) from ch_protein_coding_genes
  
//   output:
//   file("*")
  
//   script:
//   """
//   refine_orf.py \
//   --orfs $best_orfs \
//   --pb_fasta $sample_fasta \
//   --redundant ${params.sample}_redundant_accessions.txt \
//   --combined_tsv ${params.sample}_orf_combined.tsv \
//   --combined_fasta ${params.sample}_orf_combined.fasta \
//   --agg_tsv ${params.sample}_orf_aggregated.tsv \
//   --agg_fasta ${params.sample}_orf_aggregated.fasta \
//   --protein_coding_only ${params.protein_coding_only} \
//   --protein_coding_genes $ch_protein_coding_genes \
//   --cutoff ${params.refine_cutoff} \
  
//   """
// }
/*--------------------------------------------------
PacBio CDS GTF 
---------------------------------------------------*/


/*--------------------------------------------------
MetaMorpheus
---------------------------------------------------*/


/*--------------------------------------------------
Protein Inference Analysis
---------------------------------------------------*/

/*--------------------------------------------------
Peptide Analysis
---------------------------------------------------*/









