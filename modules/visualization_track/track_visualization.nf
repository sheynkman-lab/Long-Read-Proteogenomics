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
    .set { ch_sample_gtf_cds }   

Channel
    .value(file(params.refined_info))
    .ifEmpty { error "Cannot find gtf file for parameter --ch_refined_info: ${params.refined_info}" }
    .set { ch_refined_info_cds }  

Channel
    .value(file(params.refined_fasta))
    .ifEmpty { error "Cannot find gtf file for parameter --refined_db_fasta: ${params.refined_fasta}" }
    .set { ch_refined_fasta_peptide_gtf }  
    

Channel
    .value(file(params.best_orf))
    .ifEmpty { error "Cannot find gtf file for parameter --orf_calls: ${params.best_orf}" }
    .set { ch_best_orf_cds }   

Channel
    .value(file(params.pb_gene))
    .ifEmpty { error "Cannot find gtf file for parameter --pb_gene: ${params.pb_gene}" }
    .set { ch_pb_gene }   

ch_pb_gene.into{
  ch_pb_gene_cds
  ch_pb_gene_peptide_gtf
}

Channel
    .value(file(params.peptides))
    .ifEmpty { error "Cannot find gtf file for parameter --peptides: ${params.peptides}" }
    .set { ch_pacbio_peptides }  

  Channel
  .value(file(params.gencode_gtf))
    .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.gencode_gtf}" }
    .set { ch_gencode_gtf } 



process make_pacbio_cds_gtf {
  cpus 1

  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
    file(sample_gtf) from ch_sample_gtf_cds
    file(refined_info) from ch_refined_info_cds
    file(called_orfs) from ch_best_orf_cds
    file(pb_gene) from ch_pb_gene_cds
  
  output:
    file("${params.name}_with_cds.gtf") into ch_pb_cds
    file("*")
  
  script:
  """
  make_pacbio_cds_gtf.py \
  --name ${params.name} \
  --sample_gtf $sample_gtf \
  --refined_database $refined_info \
  --called_orfs $called_orfs \
  --pb_gene $pb_gene \
  --include_transcript yes

  make_pacbio_cds_gtf.py \
  --name ${params.name}_no_transcript \
  --sample_gtf $sample_gtf \
  --refined_database $refined_info \
  --called_orfs $called_orfs \
  --pb_gene $pb_gene \
  --include_transcript no
  """
}
ch_pb_cds.into{
  ch_pb_cds_bed
  ch_pb_cds_multiregion
  ch_pb_cds_peptide_gtf
  
}

/*--------------------------------------------------
Convert PacBio CDS to Bed12
---------------------------------------------------*/
process pb_cds_to_bed12 {
  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
    file(pb_cds) from ch_pb_cds_bed
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
    file(sample_gtf) from ch_pb_cds_multiregion
    file(reference_gtf) from ch_gencode_gtf
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

  when:
    params.mass_spec != false


  input:
    file(sample_gtf) from ch_pb_cds_peptide_gtf
    file(peptides) from ch_pacbio_peptides
    file(pb_gene) from ch_pb_gene_peptide_gtf
    file(refined_fasta) from ch_refined_fasta_peptide_gtf
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

  when:
    params.mass_spec != false

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











