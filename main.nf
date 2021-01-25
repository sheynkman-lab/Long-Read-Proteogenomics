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






Channel
    .value(file(params.gencode_gtf))
    .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.gencode_gtf}" }
    .set { ch_gencode_gtf }   
    
Channel
    .value(file(params.gencode_transcript_fasta))
    .ifEmpty { error "Cannot find any gencode_fasta file for parameter --gencode_fasta: ${params.gencode_transcript_fasta}" }
    .set { ch_gencode_transcript_fasta }  

Channel
    .value(file(params.gencode_translation_fasta))
    .ifEmpty { error "Cannot find any gencode_fasta file for parameter --gencode_fasta: ${params.gencode_translation_fasta}" }
    .set { ch_gencode_translation_fasta }  

Channel
    .value(file(params.sample_ccs))
    .ifEmpty { error "Cannot find file for parameter --sample_ccs: ${params.sample_ccs}" }
    .set { ch_sample_ccs }   
    
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

Channel
    .value(file(params.sample_kallisto_tpm))
    .ifEmpty { error "Cannot find any logit model file for parameter --sample_kallisto_tpm: ${params.sample_kallisto_tpm}" }
    .set { ch_sample_kallisto } 
  
Channel
    .value(file(params.normalized_ribo_kallisto))
    .ifEmpty { error "Cannot find any logit model file for parameter --normalized_ribo_kallisto: ${params.normalized_ribo_kallisto}" }
    .set { ch_normalized_ribo_kallisto } 


/*--------------------------------------------------
Reference Tables 
---------------------------------------------------*/
process generate_reference_tables {
  tag "${gencode_gtf}, ${gencode_transcript_fasta}"
  cpus 1
  publishDir "${params.outdir}/reference_tables/", mode: 'copy'

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
process make_gencode_database {
  tag "${gencode_translation_fasta}"
  cpus 1
  publishDir "${params.outdir}/gencode_db/", mode: 'copy'

  input:
  file(gencode_translation_fasta) from ch_gencode_translation_fasta
  
  output:
  file("gencode.fasta") into ch_gencode_fasta_single
  file("gencode_isoname_clusters.tsv") into ch_gencode_isoname_clusters
  
  script:
  """
  make_gencode_database.py \
  --gencode_fasta $gencode_translation_fasta \
  --output_fasta gencode.fasta \
  --output_cluster gencode_isoname_clusters.tsv
  """
}


/*--------------------------------------------------
IsoSeq3
---------------------------------------------------*/

process isoseq3 {
  tag "${sample_ccs}, ${gencode_fasta}, ${primers_fasta}"
  cpus params.max_cpus
  publishDir "${params.outdir}/isoseq3/", mode: 'copy'

  input:
  file(sample_ccs) from ch_sample_ccs
  file(gencode_fasta) from ch_gencode_fasta
  file(primers_fasta) from ch_primers_fasta
  
  output:
  file("*")
  file("${params.name}.collapsed.gff") into ch_isoseq_gtf
  // file("${params.name}.collapsed.fasta") into ch_isoseq_fasta
  file("${params.name}.collapsed.abundance.txt") into ch_fl_count
  script:
  """
  # create an index
  pbindex $sample_ccs

  # module load isoseqenv
  lima --isoseq --dump-clips --peek-guess -j ${task.cpus} $sample_ccs $primers_fasta ${params.name}.demult.bam
  isoseq3 refine --require-polya ${params.name}.demult.NEB_5p--NEB_3p.subreadset.xml $primers_fasta ${params.name}.flnc.bam

  # clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
  isoseq3 cluster ${params.name}.flnc.bam ${params.name}.polished.bam --verbose --use-qvs

  # align reads to the genome, takes few minutes (40 core machine)
  pbmm2 align $gencode_fasta ${params.name}.polished.transcriptset.xml ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${task.cpus} --log-level INFO

  # collapse redundant reads
  isoseq3 collapse ${params.name}.aligned.bam ${params.name}.collapsed.gff
  """
}
/*
Channel
  .value(file(params.fl_count))
  .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.fl_count}" }
  .set { ch_fl_count }  

Channel
  .value(file(params.sample_gtf))
  .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.sample_gtf}" }
  .set { ch_sample_gtf } 
Channel
  .value(file(params.sample_fasta))
  .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.sample_fasta}" }
  .set { ch_sample_fasta } 
*/

/*--------------------------------------------------
SQANTI3
---------------------------------------------------*/

process sqanti3 {
  tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"
  cpus params.max_cpus
  publishDir "${params.outdir}/sqanti3/", mode: 'copy'

  input:
  
  file(fl_count) from ch_fl_count
  file(gencode_gtf) from ch_gencode_gtf
  file(gencode_fasta) from ch_gencode_fasta
  file(sample_gtf) from ch_isoseq_gtf
  
  
  output:
  file("${params.name}_classification.txt") into ch_sample_classification
  file("${params.name}_corrected.fasta") into ch_sample_fasta
  file("${params.name}_corrected.gtf") into ch_sample_gtf
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

/*
Channel
  .value(file(params.sample_classification))
  .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.sample_classification}" }
  .set { ch_sample_classification } 

 */ 


/*--------------------------------------------------
Six-Frame Translation
---------------------------------------------------*/
process make_pacbio_6frm_gene_grouped {
    cpus 1
    tag "${classification}, ${ensg_gene}"
    publishDir "${params.outdir}/pacbio_6frm_gene_grouped/", mode: 'copy'

    input:
    file(classification) from ch_sample_classification
    file(ensg_gene) from ch_ensg_gene
    file(sample_fasta) from ch_sample_fasta

    output:
    file("${params.name}.6frame.fasta") into ch_6frm

    script:
    """
    make_pacbio6frm_gene_grouped.py \
    --iso_annot $classification \
    --ensg_gene $ensg_gene \
    --sample_fasta $sample_fasta \
    --output_fasta ${params.name}.6frame.fasta
    """
}


/*--------------------------------------------------
Transcriptome Summary 
---------------------------------------------------*/
process transcriptome_summary {
  cpus 1
  publishDir "${params.outdir}/transcriptome_summary/", mode: 'copy'

  input:
  file(sqanti_classification) from ch_sample_classification
  file(tpm) from ch_sample_kallisto
  file(ribo) from ch_normalized_ribo_kallisto
  file(ensg_to_gene) from ch_ensg_gene
  file(enst_to_isoname) from ch_enst_isoname
  file(len_stats) from ch_gene_lens
  
  
  output:
  file("gene_level_tab.tsv") into ch_gene_level
  file("sqanti_isoform_info.tsv") into ch_sqanti_isoform_info
  file("pb_gene.tsv") into ch_pb_gene
  
  script:
  """
  transcriptome_summary.py \
  --sq_out $sqanti_classification \
  --tpm $tpm \
  --ribo $ribo \
  --ensg_to_gene $ensg_to_gene \
  --enst_to_isoname $enst_to_isoname \
  --len_stats $len_stats 
  """
}

/*--------------------------------------------------
CPAT
---------------------------------------------------*/
process cpat {
  cpus 1
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
process orf_calling {
  tag "${orf_coord}, ${gencode_gtf}, ${sample_gtf}, ${pb_gene}, ${classification}, ${sample_fasta} "
  cpus params.max_cpus
  publishDir "${params.outdir}/orf_calling/", mode: 'copy'

  input:
  file(cpat_orfs) from ch_cpat_orfs
  file(gencode_gtf) from ch_gencode_gtf
  file(sample_gtf) from ch_sample_gtf
  file(sample_fasta) from ch_sample_fasta
  file(pb_gene) from ch_pb_gene
  file(classification) from ch_sample_classification
  
  
  output:
  file("${params.name}_best_orf.tsv") into ch_best_orf
  
  script:
  """
  orf_calling.py \
  --orf_coord $cpat_orfs \
  --gencode $gencode_gtf \
  --sample_gtf $sample_gtf \
  --pb_gene $pb_gene \
  --classification $classification \
  --sample_fasta $sample_fasta \
  --num_cores ${task.cpus} \
  --output ${params.name}_best_orf.tsv 
  """
}

/*--------------------------------------------------
Refined DB Generation 
---------------------------------------------------*/
process generate_refined_database {
  cpus 1
  tag "${best_orfs}, ${sample_fasta}, ${params.protein_coding_only}, ${protein_coding_genes}, ${params.refine_cutoff}" 

  publishDir "${params.outdir}/refined_database/", mode: 'copy'

  input:
  file(best_orfs) from ch_best_orf
  file(sample_fasta) from ch_sample_fasta
  file(protein_coding_genes) from ch_protein_coding_genes
  
  output:
  file("*")
  file("${params.name}_orf_aggregated.tsv") into ch_agg_orfs
  file("${params.name}_orf_aggregated.fasta") into ch_agg_fasta
  
  script:
  """
  refine_orf.py \
  --orfs $best_orfs \
  --pb_fasta $sample_fasta \
  --redundant ${params.name}_redundant_accessions.txt \
  --combined_tsv ${params.name}_orf_combined.tsv \
  --combined_fasta ${params.name}_orf_combined.fasta \
  --agg_tsv ${params.name}_orf_aggregated.tsv \
  --agg_fasta ${params.name}_orf_aggregated.fasta \
  --protein_coding_only ${params.protein_coding_only} \
  --protein_coding_genes $protein_coding_genes \
  --coding_score_cutoff ${params.refine_cutoff} \
  
  """
}

/*--------------------------------------------------
PacBio CDS GTF 
---------------------------------------------------*/

process make_pacbio_cds_gtf {
  cpus 1

  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
  file(sample_gtf) from ch_sample_gtf
  file(agg_orfs) from ch_agg_orfs
  file(refined_orfs) from ch_best_orf
  file(pb_gene) from ch_pb_gene
  
  output:
  file("${params.name}_cds.gtf") into ch_orf_cds
  
  script:
  """
  make_pacbio_cds_gtf.py \
  --sample_gtf $sample_gtf \
  --agg_orfs $agg_orfs \
  --refined_orfs $refined_orfs \
  --pb_gene $pb_gene \
  --output_cds ${params.name}_cds.gtf
  """
}


/*--------------------------------------------------
MetaMorpheus
---------------------------------------------------*/

// process metamorpheus {
//   publishDir "${params.outdir}/metamorpheus"

//   // input:

//   output:
//   file("*")

//   script: 
//   """
//   dotnet /metamorpheus/CMD.dll --test -v minimal -o metamorpheus
//   """
// }


/*--------------------------------------------------
Accession Mapping 
---------------------------------------------------*/

/*--------------------------------------------------
Protein Inference Analysis
---------------------------------------------------*/

/*--------------------------------------------------
Peptide Analysis
---------------------------------------------------*/









