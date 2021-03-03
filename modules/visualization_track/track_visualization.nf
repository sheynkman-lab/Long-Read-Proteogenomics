#!/usr/bin/env nextflow
/*
 * Copyright (c) 2020, Sheynkman Lab and the authors.
 *
 * This is part of the Long Read Proteogenomics Pipeline
 *
 * @authors
 * Gloria Sheynkman
 * Ben Jordan
 */

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

log.info "LRPG Visualization - N F  ~  version 0.1"
log.info "====================================="


Channel
    .value(file(params.sample_gtf))
    .ifEmpty { error "Cannot find gtf file for parameter --sample_gtf: ${params.sample_gtf}" }
    .set { ch_sample_gtf }   

Channel
    .value(file(params.refined_db_data))
    .ifEmpty { error "Cannot find gtf file for parameter --refined_db_data: ${params.refined_db_data}" }
    .set { ch_refined_db_data }  

Channel
    .value(file(params.refined_db_fasta))
    .ifEmpty { error "Cannot find gtf file for parameter --refined_db_fasta: ${params.refined_db_fasta}" }
    .set { ch_refined_db_fasta }  

Channel
    .value(file(params.orf_calls))
    .ifEmpty { error "Cannot find gtf file for parameter --orf_calls: ${params.orf_calls}" }
    .set { ch_orf_calls }   

Channel
    .value(file(params.pb_gene))
    .ifEmpty { error "Cannot find gtf file for parameter --pb_gene: ${params.pb_gene}" }
    .set { ch_pb_gene }   

Channel
    .value(file(params.peptides))
    .ifEmpty { error "Cannot find gtf file for parameter --refined_db_data: ${params.peptides}" }
    .set { ch_peptides }  

  Channel
  .value(file(params.reference_gtf))
    .ifEmpty { error "Cannot find gtf file for parameter --reference_gtf: ${params.reference_gtf}" }
    .set { ch_reference_gtf } 


/*--------------------------------------------------
PacBio CDS GTF 
---------------------------------------------------*/

process make_pacbio_cds_gtf {
  cpus 1

  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
    file(sample_gtf) from ch_sample_gtf
    file(agg_orfs) from ch_refined_db_data
    file(refined_orfs) from ch_orf_calls
    file(pb_gene) from ch_pb_gene
  
  output:
    file("${params.name}_with_cds.gtf") into ch_pb_cds
  
  script:
  """
  make_pacbio_cds_gtf.py \
  --name ${params.name} \
  --sample_gtf $sample_gtf \
  --agg_orfs $agg_orfs \
  --refined_orfs $refined_orfs \
  --pb_gene $pb_gene \
  --include_transcript ${params.include_transcript}
  """
}

/*--------------------------------------------------
Convert PacBio CDS to Bed12
---------------------------------------------------*/
process pb_cds_to_bed12 {
  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
    file(pb_cds) from ch_pb_cds
  output:
    file("${params.name}_with_cds.bed12") into ch_cds_bed
  
  script:
    """
    gtfToGenePred $pb_cds ${params.name}_with_cds.genePred
    genePredToBed ${params.name}_with_cds.genePred ${params.name}_with_cds.bed12
    """
}

/*--------------------------------------------------
Add RGB shading to tracks
---------------------------------------------------*/
process add_shading_to_cds{
  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
    file(cds_bed) from ch_cds_bed
  output:
    file("${params.name}_cds_shaded.bed12") into ch_cds_shaded

  script:
  """
  add_rgb_shading_to_pb_track.py \
  --name ${params.name} \
  --bed_file $cds_bed
  """
}

/*--------------------------------------------------
Make Region Bed for UCSC Browser
---------------------------------------------------*/
process make_multiregion{
  input:
    file(sample_gtf) from ch_pb_cds
    file(reference_gtf) from ch_reference_gtf
  output:
    file("*")

  script:
  """
  make_region_bed_for_ucsc.py \
  --name ${params.name} \
  --sample_gtf $sample_gtf \
  --reference_gtf $reference_gtf
  """
}

/*--------------------------------------------------
Make Peptide GTF 
---------------------------------------------------*/
process make_peptide_gtf{
  publishDir "${params.outdir}/peptide_track/", mode: 'copy'

  input:
    file(sample_gtf) from ch_pb_cds
    file(peptides) from ch_peptides
    file(pb_gene) from ch_pb_gene
    file(refined_fasta) from ch_refined_db_fasta
  output:
    file("${params.name}_peptides.gtf") into ch_peptide_gtf

  script:
  """
  make_peptide_gtf_file.py \
  --name ${params.name} \
  --sample_gtf $sample_gtf \
  --peptides $peptides \
  --pb_gene $pb_gene \
  --refined_fasta $refined_fasta
  """
}

/*--------------------------------------------------
Convert Peptide GTF to BED and Add RGB
---------------------------------------------------*/
process peptide_gtf_to_bed{
  publishDir "${params.outdir}/peptide_track/", mode: 'copy'

  input:
    file(peptide_gtf) from ch_peptide_gtf
    
  output:
    file("${params.name}_peptides.bed12") into ch_peptide_bed

  script:
  """
  gtfToGenePred $peptide_gtf ${params.name}_peptides.genePred
  genePredToBed ${params.name}_peptides.genePred ${params.name}_peptides.bed12

  # add rgb to colorize specific peptides 
  echo 'track name=peptide_w_specificity itemRgb=On' > ${params.name}_peptide_rgb.bed12
  cat ${params.name}_peptides.bed12 | awk '{ if (\$4 ~ /-1\$/) {\$9="0,102,0"; print \$0} else {\$9="0,51,0"; print \$0} }' >> ${params.name}_peptide_rgb.bed12
  """
  
}











