#!/usr/bin/env nextflow
/*
 * Copyright (c) 2020, Sheynkman Lab and the authors.
 *
 *   This file is part of 'proteogenomics-nf' a pipeline repository to run
 *   the "long read proteogenomics" pipeline.
 *
 *
 * @authors
 * Ben Jordan
 * Rachel Miller
 * Gloria Sheynkman
 * Christina Chatzipantsiou
 * Anne Deslattes Mays
 */

def helpMessage() {
    log.info logHeader()
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
// Header log info
log.info "\nPARAMETERS SUMMARY"
log.info "coding_score_cutoff                   : ${params.coding_score_cutoff}"
log.info "config                                : ${params.config}"
log.info "fastq_read_1                          : ${params.fastq_read_1}"
log.info "fastq_read_2                          : ${params.fastq_read_2}"
log.info "gencode_gtf                           : ${params.gencode_gtf}"
log.info "gencode_transcript_fasta              : ${params.gencode_transcript_fasta}"
log.info "gencode_translation_fasta             : ${params.gencode_translation_fasta}"
log.info "genome_fasta                          : ${params.genome_fasta}"
log.info "hexamer                               : ${params.hexamer}"
log.info "logit_model                           : ${params.logit_model}"
log.info "mainScript                            : ${params.mainScript}"
log.info "mass_spec                             : ${params.mass_spec}"
log.info "max_cpus                              : ${params.max_cpus}"
log.info "metamorpheus toml                     : ${params.metamorpheus_toml}"
log.info "name                                  : ${params.name}"
log.info "normalized_ribo_kallisto              : ${params.normalized_ribo_kallisto}"
log.info "outdir                                : ${params.outdir}"
log.info "primers_fasta                         : ${params.primers_fasta}"
log.info "rescue and resolve toml               : ${params.rescue_resolve_toml}"
log.info "sample_ccs                            : ${params.sample_ccs}"
log.info "sample_kallisto_tpm                   : ${params.sample_kallisto_tpm}"
log.info "sqanti classification                 : ${params.sqanti_classification}"
log.info "sqanti fasta                          : ${params.sqanti_fasta}"
log.info "sqanti gtf                            : ${params.sqanti_gtf}"
log.info "star_genome_dir                       : ${params.star_genome_dir}"
log.info "uniprot_protein_fasta                 : ${params.uniprot_protein_fasta}"
log.info ""


if (!params.gencode_gtf) exit 1, "Cannot find gtf file for parameter --gencode_gtf: ${params.gencode_gtf}"
ch_gencode_gtf = Channel.value(file(params.gencode_gtf))

if (!params.gencode_transcript_fasta) exit 1, "Cannot find any file for parameter --gencode_transcript_fasta: ${params.gencode_transcript_fasta}"

if (!params.gencode_translation_fasta) exit 1, "Cannot find any file for parameter --gencode_translation_fasta: ${params.gencode_translation_fasta}"

if (!params.uniprot_protein_fasta) exit 1, "Cannot find any file for parameter --uniprot_protein_fasta: ${params.uniprot_protein_fasta}"

if (!params.hexamer) exit 1, "Cannot find headmer file for parameter --hexamer: ${params.hexamer}"
ch_hexamer = Channel.value(file(params.hexamer))

if (!params.logit_model) exit 1, "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}"

if (!params.sample_kallisto_tpm) exit 1, "Cannot find any sample_kallisto_tpm file for parameter --sample_kallisto_tpm: ${params.sample_kallisto_tpm}"
ch_sample_kallisto = Channel.value(file(params.sample_kallisto_tpm))

if (!params.normalized_ribo_kallisto) exit 1, "Cannot find any normalized_ribo_kallisto file for parameter --normalized_ribo_kallisto: ${params.normalized_ribo_kallisto}"
ch_normalized_ribo_kallisto = Channel.value(file(params.normalized_ribo_kallisto))

//if (!params.metamorpheus_toml) exit 1, "Cannot find any file for parameter --metamorpheus_toml: ${params.metamorpheus_toml}"

// if (!params.fastq_read_1) exit 1, "No file found for the parameter --fastq_read_1 at the location ${params.fastq_read_1}"
// if (!params.fastq_read_2) exit 1, "No file found for the parameter --fastq_read_2 at the location ${params.fastq_read_2}"
if (params.mass_spec != false & params.rescue_resolve_toml == false){
   exit 1, "Cannot find file for parameter --rescue_resolve_toml: ${params.rescue_resolve_toml}"
} else if (params.mass_spec != false & params.rescue_resolve_toml != false){
   ch_rr_toml = Channel.value(file(params.rescue_resolve_toml))
} else{
   ch_rr_toml = Channel.from("no mass spec")
}

if (params.gencode_translation_fasta.endsWith('.gz')){
   ch_gencode_translation_fasta = Channel.value(file(params.gencode_translation_fasta))
} else {
   ch_gencode_translation_fasta_uncompressed = Channel.value(file(params.gencode_translation_fasta))
}

if (params.gencode_transcript_fasta.endsWith('.gz')){
   ch_gencode_transcript_fasta = Channel.value(file(params.gencode_transcript_fasta))
} else {
   ch_gencode_transcript_fasta_uncompressed = Channel.value(file(params.gencode_transcript_fasta))
}

if (params.genome_fasta.endsWith('.gz')){
   ch_genome_fasta = Channel.value(file(params.genome_fasta))
} else {
   ch_genome_fasta_uncompressed = Channel.value(file(params.genome_fasta))
}

if (params.uniprot_protein_fasta.endsWith('.gz')){
   ch_uniprot_protein_fasta = Channel.value(file(params.uniprot_protein_fasta))
} else {
   ch_uniprot_protein_fasta_uncompressed = Channel.value(file(params.uniprot_protein_fasta))
}

if (params.logit_model.endsWith('.gz')) {
   ch_logit_model = Channel.value(file(params.logit_model))
} else {
   ch_logit_model_uncompressed = Channel.value(file(params.logit_model))
}

if (!params.sqanti_fasta == false) {
   if (params.sqanti_fasta.endsWith('.gz')) {
      ch_sqanti_fasta = Channel.value(file(params.sqanti_fasta))
   } else {
      ch_sqanti_fasta_uncompressed = Channel.value(file(params.sqanti_fasta))
   }
}

ch_fastq_reads = Channel.from(params.fastq_read_1, params.fastq_read_2).filter(String).flatMap{ files(it) }

// Implements logic for cloud compatibility, NO_TOML_FILE as variable only works for envs with local file system
projectDir = workflow.projectDir

if (!params.metamorpheus_toml) {
    log.warn "No toml file specified via --metamorpheus_toml for Metamorpheus, proceeding without"
    ch_metamorpheus_toml = Channel.value(file("${projectDir}/assets/NO_TOML_FILE"))
}

if (params.metamorpheus_toml) {
    log.warn "Metamorpheus toml file specified: ${params.metamorpheus_toml}"
    ch_metamorpheus_toml = Channel.value(file(params.metamorpheus_toml))
}

ch_metamorpheus_toml.into{
   ch_metamorpheus_toml_gencode
   ch_metamorpheus_toml_uniprot
   ch_metamorpheus_toml_pacbio_refined
   ch_metamorpheus_toml_pacbio_filtered
   ch_metamorpheus_toml_pacbio_hybrid
}

if (!params.mass_spec == false) {
   if (!params.mass_spec.endsWith("tar.gz")) {
      ch_mass_spec_raw = Channel.fromPath("${params.mass_spec}/*.raw")
      ch_mass_spec_mzml = Channel.fromPath("${params.mass_spec}/*.{mzml,mzML}")
   } else {
      if (params.mass_spec.endsWith("tar.gz")){
         ch_mass_spec_raw_mzml_tar_gz = Channel.value(file(params.mass_spec))
      }
   }
} else {
   ch_mass_spec_raw = Channel.from("no mass spec")
   ch_mass_spec_mzml = Channel.from("no mass spec")
}

if (!params.star_genome_dir == false) {
   if (!params.star_genome_dir.endsWith("tar.gz")) {
      ch_genome_dir = Channel.fromPath(params.star_genome_dir, type:'dir')
   } else {
      if (params.star_genome_dir.endsWith("tar.gz")) {
         ch_genome_dir_tar_gz = Channel.fromPath(params.star_genome_dir)
      }
   }
}

/*--------------------------------------------------
Decompress Logit Model
---------------------------------------------------*/
if (params.logit_model.endsWith('.gz')) {
   process gunzip_logit_model {
      tag "decompress logit model"
      cpus 1

      input:
      file(logit_model) from ch_logit_model

      output:
      file("*.RData") into ch_logit_model_uncompressed

      script:
      """
      gunzip -f ${logit_model}
      """
   }
}

/*--------------------------------------------------
Untar & decompress mass spec file
---------------------------------------------------*/
if (params.mass_spec != false) {
   if (params.mass_spec.endsWith("tar.gz")) {
      process untar_mass_spec {
         tag "${raw_mzml_tar_gz}"
         cpus 1

         input:
         file(raw_mzml_tar_gz) from ch_mass_spec_raw_mzml_tar_gz

         output:
         file("${raw_mzml_tar_gz.simpleName}/*.raw") optional true into ch_mass_spec_raw
         file("${raw_mzml_tar_gz.simpleName}/*.{mzml,mzML}") optional true into ch_mass_spec_mzml

         script:
         """
         tar xvzf $raw_mzml_tar_gz
         """
      }
   }
}

/*--------------------------------------------------
Decompress gencode translation fasta file
---------------------------------------------------*/
if (params.gencode_translation_fasta.endsWith('.gz')) {
   process gunzip_gencode_translation_fasta {
   tag "decompress gzipped gencode translation fasta"
   cpus 1

   input:
   file(gencode_translation_fasta) from ch_gencode_translation_fasta

   output:
   file("*.{fa,fasta}") into ch_gencode_translation_fasta_uncompressed

   script:
   """
   gunzip -f ${gencode_translation_fasta}
   """
   }
}

/*--------------------------------------------------
Decompress gencode transcript fasta file
---------------------------------------------------*/
if (params.gencode_transcript_fasta.endsWith('.gz')) {
   process gunzip_gencode_transcript_fasta {
   tag "decompress gzipped gencode transcript fasta"
   cpus 1

   input:
   file(gencode_transcript_fasta) from ch_gencode_transcript_fasta

   output:
   file("*.{fa,fasta}") into ch_gencode_transcript_fasta_uncompressed

   script:
   """
   gunzip -f ${gencode_transcript_fasta}
   """
   }
}

/*--------------------------------------------------
Decompress genome fasta file
---------------------------------------------------*/
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

/*--------------------------------------------------
Decompress uniprot protein fasta file
---------------------------------------------------*/
if (params.uniprot_protein_fasta.endsWith('.gz')) {
   process gunzip_uniprot_protein_fasta {
   tag "decompress gzipped uniprot protein fasta"
   cpus 1

   input:
   file(uniprot_protein_fasta) from ch_uniprot_protein_fasta

   output:
   file("*.{fa,fasta}") into ch_uniprot_protein_fasta_uncompressed

   script:
   """
   gunzip -f ${uniprot_protein_fasta}
   """
   }
}

/*--------------------------------------------------
Reference Tables 
---------------------------------------------------*/

process generate_reference_tables {
  tag "${gencode_gtf}, ${gencode_transcript_fasta}"
  cpus 1
  publishDir "${params.outdir}/${params.name}/reference_tables/", mode: 'copy'

  input:
  file(gencode_gtf) from ch_gencode_gtf
  file(gencode_transcript_fasta) from ch_gencode_transcript_fasta_uncompressed
  
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

// partition channels for use by multiple modules
ch_protein_coding_genes.into{
  ch_protein_coding_genes_gencode_fasta
  ch_protein_coding_genes_filter_sqanti
}

ch_ensg_gene.into{
  ch_ensg_gene_filter
  ch_ensg_gene_six_frame
  ch_ensg_gene_pclass
}
ch_gene_lens.into{
  ch_gene_lens_transcriptome
  ch_gene_lens_aggregate
}

ch_gene_isoname.into{
  ch_gene_isoname_pep_viz
  ch_gene_isoname_pep_analysis
}



/*--------------------------------------------------
GENCODE Database
 * Clusters same-protein GENCODE entries 
 * Output:  gencode_protein.fasta
 *            - gencode fasta file without duplicate entries
 *          gencode_isoname_clusters.tsv
 *            - transcript accessions that were clustered
---------------------------------------------------*/
process make_gencode_database {
  tag "${gencode_translation_fasta}"
  cpus 1
  publishDir "${params.outdir}/${params.name}/gencode_db/", mode: 'copy'

  input:
  file(gencode_translation_fasta) from ch_gencode_translation_fasta_uncompressed
  file(protein_coding_genes) from ch_protein_coding_genes_gencode_fasta
  output:
  file("gencode_protein.fasta") into ch_gencode_protein_fasta
  file("gencode_isoname_clusters.tsv")
  
  script:
  """
  make_gencode_database.py \
  --gencode_fasta $gencode_translation_fasta \
  --protein_coding_genes $protein_coding_genes \
  --output_fasta gencode_protein.fasta \
  --output_cluster gencode_isoname_clusters.tsv
  """
}

ch_gencode_protein_fasta.into{
  ch_gencode_protein_fasta_metamorpheus
  ch_gencode_protein_fasta_mapping
  ch_gencode_protein_fasta_hybrid
  ch_gencode_protein_fasta_novel

}

/*--------------------------------------------------
Run IsoSeq3 and SQANTI3 if SQANTI output not provided
---------------------------------------------------*/
if( params.sqanti_classification==false || params.sqanti_fasta==false || params.sqanti_gtf==false ){
  if (!params.sample_ccs) exit 1, "Cannot find file for parameter --sample_ccs: ${params.sample_ccs}"
  ch_sample_ccs = Channel.value(file(params.sample_ccs))

  if (!params.primers_fasta) exit 1, "Cannot find any seq file for parameter --primers_fasta: ${params.primers_fasta}"
  ch_primers_fasta = Channel.value(file(params.primers_fasta))

  if (!params.genome_fasta) exit 1, "Cannot find any seq file for parameter --genome_fasta: ${params.genome_fasta}"
  
  ch_genome_fasta_uncompressed.into{
    ch_genome_fasta_star
    ch_genome_fasta_isoseq
    ch_genome_fasta_sqanti
  }
  
  /*--------------------------------------------------
  IsoSeq3
   * Runs IsoSeq3 on CCS reads, aligning to genome and
   * collapsing redundant reads
   * STEPS
   *   - ensure only qv10 reads from ccs are kept as input
   *   - find and remove adapters/barcodes
   *   - filter for non-concatamer, polya containing reads
   *   - clustering of reads
   *   - align reads to the genome
   *   - collapse redundant reads
  ---------------------------------------------------*/
  process isoseq3 {
    tag "${sample_ccs}, ${genome_fasta}, ${primers_fasta}"
    cpus params.max_cpus
    publishDir "${params.outdir}/${params.name}/isoseq3/", mode: 'copy'

    input:
    file(sample_ccs) from ch_sample_ccs
    file(genome_fasta) from ch_genome_fasta_isoseq
    file(primers_fasta) from ch_primers_fasta
    
    output:
    file("${params.name}.collapsed.gff") into ch_isoseq_gtf
    file("${params.name}.collapsed.abundance.txt") into ch_fl_count
    file("${params.name}.collapsed.fasta")
    file("${params.name}.collapsed.report.json")
    file("${params.name}.demult.lima.summary")
    file("${params.name}.flnc.bam")
    file("${params.name}.flnc.bam.pbi")
    file("${params.name}.flnc.filter_summary.json")


    script:
    """
    # ensure that only qv10 reads from ccs are input
    bamtools filter -tag 'rq':'>=0.90' -in $sample_ccs -out filtered.$sample_ccs 

    # create an index for the ccs bam
    pbindex filtered.$sample_ccs

    # find and remove adapters/barcodes
    lima --isoseq --dump-clips --peek-guess -j ${task.cpus} filtered.$sample_ccs $primers_fasta ${params.name}.demult.bam

    # filter for non-concatamer, polya containing reads
    isoseq3 refine --require-polya ${params.name}.demult.NEB_5p--NEB_3p.bam $primers_fasta ${params.name}.flnc.bam

    # clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
    isoseq3 cluster ${params.name}.flnc.bam ${params.name}.clustered.bam --verbose --use-qvs

    # align reads to the genome, takes few minutes (40 core machine)
    pbmm2 align $genome_fasta ${params.name}.clustered.hq.bam ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${task.cpus} --log-level INFO

    # collapse redundant reads
    isoseq3 collapse ${params.name}.aligned.bam ${params.name}.collapsed.gff
    """
  }

  /*--------------------------------------------------
   Untar & Decompress star genome directory
  ---------------------------------------------------*/
  if (params.star_genome_dir != false) {
     if (params.star_genome_dir.endsWith("tar.gz")) {

         process untar_star_genome_dir {
            tag "${genome_dir_tar_gz}"
            cpus 1

            input:
            file(genome_dir_tar_gz) from ch_genome_dir_tar_gz

            output:
            file("star_genome") into ch_genome_dir

            script:
            """
            tar xvzf $genome_dir_tar_gz
            """
        }
     }
  }

  /*--------------------------------------------------
  STAR Alignment
   * STAR alignment is run only if sqanti has not been
   *  previously been run and if fastq (short read RNAseq) files
   *  have been provided.
   *  if( params.sqanti_classification==false || params.sqanti_fasta==false || params.sqanti_gtf==false )
   *  STAR alignment is run if fastq reads are provided
   *  Junction alignments are fed to SQANTI3 where 
   *  information is used in classificaiton filtering
   *
   * STEPS
   *   - generate star genome index (skipped if provided) 
   *   - star read alignment 
  ---------------------------------------------------*/
  if(!params.star_genome_dir){
    process star_generate_genome{
          cpus params.max_cpus
          publishDir "${params.outdir}/${params.name}/star_index", mode: "copy"
          
          when:
          (params.fastq_read_1 != false | params.fastq_read_2 !=false) & params.star_genome_dir == false

          input :
              file(gencode_gtf) from ch_gencode_gtf
              file(genome_fasta) from ch_genome_fasta_star

          output:
              path("star_genome") into ch_genome_dir

          script:
          """
          mkdir star_genome
          STAR --runThreadN  ${task.cpus} \
          --runMode genomeGenerate \
          --genomeDir star_genome \
          --genomeFastaFiles $genome_fasta \
          --sjdbGTFfile $gencode_gtf \
          --genomeSAindexNbases 11
          """
      }
  }

  if(params.fastq_read_1 != false | params.fastq_read_2 !=false){
      process star_alignment{
          cpus params.max_cpus
          publishDir "${params.outdir}/${params.name}/star", mode: "copy"
          when:
              params.fastq_read_1 != false | params.fastq_read_2 !=false

          input :
              file(fastq_reads) from ch_fastq_reads.collect()
              path(genome_dir) from ch_genome_dir

          output:
              file("*SJ.out.tab") into ch_star_junction
              file("*Log.final.out")

          script:
          """
          STAR --runThreadN ${task.cpus} \
          --genomeDir $genome_dir \
          --outFileNamePrefix ./${params.name} \
          --readFilesIn $fastq_reads \
          --readFilesCommand zcat
          """
      }
  }
  else{
      Channel
          .from('none')
          .set{ch_star_junction}
  }


  /*--------------------------------------------------
  SQANTI3
   * https://github.com/ConesaLab/SQANTI3
   * Corrects any errors in alignment from IsoSeq3 and 
   * classifies each accession in relation to the reference
   * genome
  ---------------------------------------------------*/
  process sqanti3 {
    tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"
    cpus params.max_cpus
    publishDir "${params.outdir}/${params.name}/sqanti3/", mode: 'copy'

    input:
    file(fl_count) from ch_fl_count
    file(gencode_gtf) from ch_gencode_gtf
    file(genome_fasta) from ch_genome_fasta_sqanti
    file(sample_gtf) from ch_isoseq_gtf
    file(star_junction) from ch_star_junction
    
    
    output:
    file("${params.name}_classification.txt") into ch_sample_unfiltered_classification
    file("${params.name}_corrected.fasta") into ch_sample_unfiltered_fasta
    file("${params.name}_corrected.gtf") into ch_sample_unfiltered_gtf
    file("${params.name}_junctions.txt")
    file("${params.name}_sqanti_report.pdf")
    file("${params.name}.params.txt")
    
    script:
    if(params.fastq_read_1 == false & params.fastq_read_2 ==false)
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
    else
      """
      sqanti3_qc.py \
      $sample_gtf \
      $gencode_gtf \
      $genome_fasta \
      --skipORF \
      -o ${params.name} \
      --fl_count $fl_count  \
      --gtf \
      -c $star_junction
      """
    //
  }


} 
else{
  ch_sample_unfiltered_classification = Channel.value(file(params.sqanti_classification))
//  ch_sample_unfiltered_fasta = Channel.value(file(params.sqanti_fasta))
  ch_sample_unfiltered_fasta = ch_sqanti_fasta_uncompressed
  ch_sample_unfiltered_gtf = Channel.value(file(params.sqanti_gtf))
}

/*--------------------------------------------------
Filter SQANTI 
 * Filter SQANTI results based on several criteria 
 * - protein coding only
 *      PB transcript aligns to a GENCODE-annotated protein coding gene.
 * - percent A downstream
 *      perc_A_downstreamTTS : percent of genomic "A"s in the downstream 20 bp window. 
 *      If this number if high (> 80%), the 3' end have arisen from intra-priming during the RT step 
 * - RTS stage
 *      RTS_stage: TRUE if one of the junctions could be an RT template switching artifact.
 * - Structural Category
 *      keep only transcripts that have a isoform structural category of:
 *        -novel_not_in_catalog
 *        -novel_in_catalog
 *        -incomplete-splice_match
 *        -full-splice_match 
---------------------------------------------------*/
process filter_sqanti {
  publishDir "${params.outdir}/${params.name}/sqanti3-filtered/", mode: 'copy'

  input:
    file(classification) from ch_sample_unfiltered_classification
    file(sample_fasta) from ch_sample_unfiltered_fasta
    file(sample_gtf) from ch_sample_unfiltered_gtf
    file(protein_coding_genes) from ch_protein_coding_genes_filter_sqanti
    file(ensg_gene) from ch_ensg_gene_filter

  output:
    file("${params.name}_classification.5degfilter.tsv") into ch_sample_classification
    file("${params.name}_corrected.5degfilter.fasta") into ch_sample_fasta
    file("${params.name}_corrected.5degfilter.gff") into ch_sample_gtf
    file("*")

  script:
    """
    filter_sqanti.py \
    --sqanti_classification $classification \
    --sqanti_corrected_fasta $sample_fasta \
    --sqanti_corrected_gtf $sample_gtf \
    --protein_coding_genes $protein_coding_genes \
    --ensg_gene $ensg_gene \
    --filter_protein_coding yes \
    --filter_intra_polyA yes \
    --filter_template_switching yes \
    --percent_A_downstream_threshold 95 \
    --structural_categories_level strict \
    --minimum_illumina_coverage 3 \

    collapse_isoforms.py \
    --name ${params.name} \
    --sqanti_gtf filtered_${params.name}_corrected.gtf \
    --sqanti_fasta filtered_${params.name}_corrected.fasta

    collapse_classification.py \
    --name ${params.name} \
    --collapsed_fasta ${params.name}_corrected.5degfilter.fasta \
    --classification filtered_${params.name}_classification.tsv
    """
}

ch_sample_classification.into{
  ch_sample_classification_six_frame
  ch_sample_classification_transcriptome
  ch_sample_classification_orf
}

ch_sample_fasta.into{
  ch_sample_fasta_cpat
  ch_sample_fasta_six_frame
  ch_sample_fasta_orf
  ch_sample_fasta_refine
}

ch_sample_gtf.into{
  ch_sample_gtf_orf
  ch_sample_gtf_cds
}


/*--------------------------------------------------
Six-Frame Translation
 * Generates a fasta file of all possible protein sequences
 * derivable from each PacBio transcript, by translating the
 * fasta file in all six frames (3+, 3-). This is used to examine
 * what peptides could theoretically match the peptides found via
 * a mass spectrometry search against GENCODE. 
---------------------------------------------------*/
process six_frame_translation {
    cpus 1
    tag "${classification}, ${ensg_gene}"
    publishDir "${params.outdir}/${params.name}/pacbio_6frm_gene_grouped/", mode: 'copy'

    input:
    file(classification) from ch_sample_classification_six_frame
    file(ensg_gene) from ch_ensg_gene_six_frame
    file(sample_fasta) from ch_sample_fasta_six_frame

    output:
    file("${params.name}.6frame.fasta") into ch_six_frame

    script:
    """
    six_frame_translation.py \
    --iso_annot $classification \
    --ensg_gene $ensg_gene \
    --sample_fasta $sample_fasta \
    --output_fasta ${params.name}.6frame.fasta
    """
}


/*--------------------------------------------------
Transcriptome Summary
 * Compares the abundance (CPM) based on long-read sequencing
 * to the abundances (TPM) inferred from short-read sequencing,
 * as computed by Kallisto (analyzed outside of this pipeline).
 * Additionally produces a pacbio-gene reference table
---------------------------------------------------*/
process transcriptome_summary {
  cpus 1
  publishDir "${params.outdir}/${params.name}/transcriptome_summary/", mode: 'copy'

  input:
  file(sqanti_classification) from ch_sample_classification_transcriptome
  file(tpm) from ch_sample_kallisto
  file(ribo) from ch_normalized_ribo_kallisto
  file(ensg_to_gene) from ch_ensg_gene
  file(enst_to_isoname) from ch_enst_isoname
  file(len_stats) from ch_gene_lens_transcriptome
  
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

ch_pb_gene.into{
  ch_pb_gene_orf
  ch_pb_gene_cds
  ch_pb_gene_peptide_analysis
  
}


/*--------------------------------------------------
CPAT
 * CPAT is a bioinformatics tool to predict an RNA’s coding probability 
 * based on the RNA sequence characteristics. 
 * To achieve this goal, CPAT calculates scores of sequence-based features 
 * from a set of known protein-coding genes and background set of non-coding genes.
 *     ORF size
 *     ORF coverage
 *     Fickett score
 *     Hexamer usage bias
 * 
 * CPAT will then builds a logistic regression model using these 4 features as 
 * predictor variables and the “protein-coding status” as the response variable. 
 * After evaluating the performance and determining the probability cutoff, 
 * the model can be used to predict new RNA sequences.
 *
 * https://cpat.readthedocs.io/en/latest/
---------------------------------------------------*/
process cpat {
  cpus 1
  tag "${hexamer}, ${logit_model}, ${sample_fasta}"

  publishDir "${params.outdir}/${params.name}/cpat/", mode: 'copy'

  input:
  file(hexamer) from ch_hexamer
  file(logit_model) from ch_logit_model_uncompressed
  file(sample_fasta) from ch_sample_fasta_cpat

  output:
  file("${params.name}.ORF_prob.tsv") into ch_cpat_all_orfs
  file("${params.name}.ORF_prob.best.tsv") into ch_cpat_best_orf
  file("${params.name}.ORF_seqs.fa") into ch_cpat_protein_fasta
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

ch_cpat_all_orfs.into{
  ch_cpat_all_orfs_for_orf_calling
  ch_cpat_all_orfs_for_peptide_analysis
}

ch_cpat_protein_fasta.into{
  ch_cpat_protein_fasta_orf_calling
  ch_cpat_protein_fasta_peptide_analysis
}

/*--------------------------------------------------
ORF Calling 
 * Selects the most plausible ORF from each pacbio transcript,
 * using the following information
 *    comparison of ATG start to reference (GENCODE) 
 *        - selects ORF with ATG start matching the one in the reference, if it exists 
 *    coding probability score from CPAT
 *    number of upstream ATGs for the candidate ORF
 *        - decrease score as number of upstream ATGs increases 
 *         using sigmoid function
 *  Additionally provides calling confidence of each ORF called
 *      - Clear Best ORF  : best score and fewest upstream ATGs of all called ORFs
 *      - Plausible ORF   : not clear best, but decent CPAT coding_score (>0.364) 
 *      - Low Quality ORF : low CPAT coding_score (<0.364)       
---------------------------------------------------*/
process orf_calling {
  tag "${orf_coord}, ${gencode_gtf}, ${sample_gtf}, ${pb_gene}, ${classification}, ${sample_fasta} "
  cpus params.max_cpus
  publishDir "${params.outdir}/${params.name}/orf_calling/", mode: 'copy'

  input:
  file(cpat_orfs) from ch_cpat_all_orfs_for_orf_calling
  file(cpat_fasta) from ch_cpat_protein_fasta_orf_calling
  file(gencode_gtf) from ch_gencode_gtf
  file(sample_gtf) from ch_sample_gtf_orf
  file(sample_fasta) from ch_sample_fasta_orf
  file(pb_gene) from ch_pb_gene_orf
  file(classification) from ch_sample_classification_orf
  
  output:
  file("${params.name}_best_orf.tsv") into ch_best_orf
  
  script:
  """
  orf_calling.py \
  --orf_coord $cpat_orfs \
  --orf_fasta $cpat_fasta \
  --gencode $gencode_gtf \
  --sample_gtf $sample_gtf \
  --pb_gene $pb_gene \
  --classification $classification \
  --sample_fasta $sample_fasta \
  --num_cores ${task.cpus} \
  --output ${params.name}_best_orf.tsv 
  """
}

ch_best_orf.into{
  ch_best_orf_refine
  ch_best_orf_cds
  ch_best_orf_sqanti_protein
  ch_best_orf_pclass
}

/*--------------------------------------------------
Refined DB Generation
 * - Filteres ORF database to only include accessions 
 *   with a CPAT coding score above a threshold (default 0.0)
 * - Filters ORFs to only include ORFs that have a stop codon 
 * - Collapses transcripts that produce the same protein
 *   into one entry, keeping a base accession (first alphanumeric).
 *   Abundances of transcripts (CPM) are collapsed during this process.
---------------------------------------------------*/
process refine_orf_database {
  cpus 1
  tag "${best_orfs}, ${sample_fasta}, ${params.coding_score_cutoff}" 

  publishDir "${params.outdir}/${params.name}/refined_database/", mode: 'copy'

  input:
  file(best_orfs) from ch_best_orf_refine
  file(sample_fasta) from ch_sample_fasta_refine
  
  
  output:
  // file("*")
  file("${params.name}_orf_refined.tsv") into ch_refined_info
  file("${params.name}_orf_refined.fasta") into ch_refined_fasta
  
  script:
  """
  refine_orf_database.py \
  --name ${params.name} \
  --orfs $best_orfs \
  --pb_fasta $sample_fasta \
  --coding_score_cutoff ${params.coding_score_cutoff} \
  """
}
ch_refined_info.into{
  ch_refined_info_rename
  ch_refined_info_cds
  ch_refined_info_pclass
  ch_refined_info_aggregate
}



/*--------------------------------------------------
PacBio CDS GTF 
 * derive a GTF file that includes the ORF regions (as CDS features)
---------------------------------------------------*/
process make_pacbio_cds_gtf {
  cpus 1

  publishDir "${params.outdir}/${params.name}/pacbio_cds/", mode: 'copy'

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
  ch_pb_cds_rename_cds
  ch_pb_cds_5p_utr
  ch_pb_cds_rename_pr
}

/*--------------------------------------------------
Rename CDS to Exon
 * Preprocessing step to SQANTI Protein
 * CDS is renamed to exon and transcript stop and start
 * locations are updated to reflect CDS start and stop
---------------------------------------------------*/
process rename_cds_to_exon{
    publishDir "${params.outdir}/${params.name}/rename_cds/", mode: 'copy'
    tag "${params.name} ${reference_gtf} ${sample_gtf}"
    cpus params.max_cpus
    input:
        file(reference_gtf) from ch_gencode_gtf
        file(sample_gtf) from ch_pb_cds_rename_cds
    output:
        // file("*")
        file("${params.name}.cds_renamed_exon.gtf") into ch_sample_cds_renamed
        file("${params.name}.transcript_exons_only.gtf") into ch_sample_transcript_exon_only
        file("gencode.cds_renamed_exon.gtf") into ch_ref_cds_renamed
        file("gencode.transcript_exons_only.gtf") into ch_ref_transcript_exon_only

    script:
        """
        rename_cds_to_exon.py \
        --sample_gtf $sample_gtf \
        --sample_name ${params.name} \
        --reference_gtf $reference_gtf \
        --reference_name gencode \
        --num_cores ${params.max_cpus}
        """
}
/*--------------------------------------------------
SQANTI Protein
 * Classify protein splice sites and calculates additional
 * statistics for start and stop of ORF
---------------------------------------------------*/
process sqanti_protein{
    publishDir "${params.outdir}/${params.name}/sqanti_protein/", mode: 'copy'
    input:
        file(sample_exon) from ch_sample_transcript_exon_only
        file(sample_cds) from ch_sample_cds_renamed
        file(reference_exon) from ch_ref_transcript_exon_only
        file(reference_cds) from ch_ref_cds_renamed
        file(best_orf) from ch_best_orf_sqanti_protein
    output:
        file("${params.name}.sqanti_protein_classification.tsv") into ch_sqanti_protein_classification
    script:
    """
    sqanti3_protein.py \
    $sample_exon \
    $sample_cds \
    $best_orf \
    $reference_exon \
    $reference_cds \
    -d ./ \
    -p ${params.name}
    """
}


/*--------------------------------------------------
5' UTR Status
 * Intermediate step for protein classification
 * Dtermines the 5' UTR status of the protein in order 
 * to classify protein category in latter step
---------------------------------------------------*/
process five_prime_utr{
  // publishDir "${params.outdir}/${params.name}/5p_utr/", mode: 'copy'
  input:
    file(reference_gtf) from ch_gencode_gtf
    file(sample_cds) from ch_pb_cds_5p_utr
    file(sqanti_protein_classification) from ch_sqanti_protein_classification
  output:
    file("${params.name}.sqanti_protein_classification_w_5utr_info.tsv") into ch_sqanti_protein_classification_w_5utr
  script:
    """
    1_get_gc_exon_and_5utr_info.py \
    --gencode_gtf $reference_gtf \
    --odir ./

    2_classify_5utr_status.py \
    --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \
    --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
    --sample_cds_gtf $sample_cds \
    --odir ./ 

    
    3_merge_5utr_info_to_pclass_table.py \
    --name ${params.name} \
    --utr_info pb_5utr_categories.tsv \
    --sqanti_protein_classification $sqanti_protein_classification \
    --odir ./
    """
}

/*--------------------------------------------------
Protein Classification
 * Classifies protein based on splicing and start site
 * main classifications are 
 *   pFSM: full-protein-match
 *     - protein fully matches a gencode protein
 *   pISM: incomplete-protein-match
 *     - protein only partially matches gencode protein
 *     - considered an N- or C-terminus truncation artifact
 *   pNIC: novel-in-catelog
 *     - protein composed of known N-term, splicing, and/or C-term in new combinations
 *   pNNC: novel-not-in-catelog
 *     - protein composed of novel N-term, splicing, and/or C-terminus
---------------------------------------------------*/
process protein_classification{
  publishDir "${params.outdir}/${params.name}/protein_classification/", mode: 'copy'
  input:
    file(protein_classification) from ch_sqanti_protein_classification_w_5utr
    file(best_orf) from ch_best_orf_pclass
    file(refined_info) from ch_refined_info_pclass
    file(ensg_gene) from ch_ensg_gene_pclass
  
  output:
    file("${params.name}_unfiltered.protein_classification.tsv") into ch_protein_classification_unfiltered
    file("${params.name}_genes.tsv") into ch_pr_genes
  script:
    """
    protein_classification_add_meta.py \
    --protein_classification $protein_classification \
    --best_orf $best_orf \
    --refined_meta $refined_info \
    --ensg_gene $ensg_gene \
    --name ${params.name} \
    --dest_dir ./


    protein_classification.py \
    --sqanti_protein ${params.name}.protein_classification_w_meta.tsv \
    --name ${params.name}_unfiltered \
    --dest_dir ./
    """
}
ch_pr_genes.into{
  ch_pr_genes_rename
  ch_pb_gene_peptide_gtf
}


/*--------------------------------------------------
Protein Gene Rename
 * Mapings of PacBio transcripts/proteins to GENCODE genes.
 * Some PacBio transcripts and the associated PacBio
 * predicted protein can map two different genes.
 * Some transcripts can also map to multiple genes.
---------------------------------------------------*/
process protein_gene_rename{
  publishDir "${params.outdir}/${params.name}/protein_gene_rename/", mode: 'copy'
  input:
    file(protein_genes) from ch_pr_genes_rename
    file(sample_cds) from ch_pb_cds_rename_pr
    file(refined_fasta) from ch_refined_fasta
    file(refined_info) from ch_refined_info_rename
  output:
    file("${params.name}_with_cds_refined.gtf") into ch_renamed_refined_cds
    file("${params.name}.protein_refined.fasta") into ch_renamed_refined_fasta
    file("${params.name}_orf_refined_gene_update.tsv") into ch_renamed_refined_info
  script:
    """
    protein_gene_rename.py \
    --sample_gtf $sample_cds \
    --sample_protein_fasta $refined_fasta \
    --sample_refined_info $refined_info \
    --pb_protein_genes $protein_genes \
    --name ${params.name}
    """ 
}
ch_renamed_refined_fasta.into{
  ch_refined_fasta_metamorpheus
  ch_refined_fasta_rescue_resolve
  ch_refined_fasta_pep_analysis
  ch_refined_fasta_peptide_gtf
  ch_refined_fasta_peptide_viz
  ch_refined_fasta_mapping
  ch_refined_fasta_pclass_filter
}

ch_renamed_refined_cds.into{
  ch_pb_cds_filter
  ch_pb_cds_bed
  ch_pb_cds_multiregion
  ch_pb_cds_peptide_gtf
}

/*--------------------------------------------------
Protein Filtering
 * Filters out proteins that are:
 *  - not pFSM, pNIC, pNNC
 *  - are pISMs (either N-terminus or C-terminus truncations)
 *  - pNNC with junctions after the stop codon (default 2)  
---------------------------------------------------*/
process filter_protein{
  publishDir "${params.outdir}/${params.name}/protein_filter/", mode: 'copy'
  input:
    file(reference_gtf) from ch_gencode_gtf
    file(protein_classification) from ch_protein_classification_unfiltered
    file(protein_fasta) from ch_refined_fasta_pclass_filter
    file(sample_cds) from ch_pb_cds_filter
  output:
   file("${params.name}.classification_filtered.tsv") into ch_filtered_protein_classification
   file("${params.name}.filtered_protein.fasta") into ch_filtered_protein_fasta
   file("${params.name}_with_cds_filtered.gtf") into ch_filtered_cds
  script:
    """
    protein_filter.py \
    --protein_classification $protein_classification \
    --gencode_gtf $reference_gtf \
    --protein_fasta $protein_fasta \
    --sample_cds_gtf $sample_cds \
    --min_junctions_after_stop_codon ${params.min_junctions_after_stop_codon} \
    --name ${params.name} \
    """
}
ch_filtered_cds.into{
  ch_filtered_cds_bed
  ch_filtered_cds_agg
  ch_filtered_cds_multiregion
}
ch_filtered_protein_fasta.into{
  ch_filtered_protein_fasta_metamorpheus
  ch_filtered_protein_fasta_aggregate
  ch_filtered_protein_fasta_pep_analysis
}

/*--------------------------------------------------
Protein Hybrid Database
 * Makes a hybrid database that is composed of 
 * high-confidence PacBio proteins and GENCODE proteins
 * for genes that are not in the high-confidence space
 * High-confidence is defined as genes in which the PacBio
 * sampling is adequate (average transcript length 1-4kb
 * and a total of 3 CPM per gene
---------------------------------------------------*/
process make_hybrid_database{
  publishDir "${params.outdir}/${params.name}/hybrid_protein_database/", mode: 'copy'
  
  input:
    file(protein_classification) from ch_filtered_protein_classification
    file(gene_lens) from ch_gene_lens_aggregate
    file(pb_fasta) from ch_filtered_protein_fasta_aggregate
    file(gc_fasta) from ch_gencode_protein_fasta_hybrid
    file(refined_info) from ch_renamed_refined_info
    file(sample_cds) from ch_filtered_cds_agg
  output:
    file("*")
    file("${params.name}_cds_high_confidence.gtf")     into ch_high_confidence_cds
    file("${params.name}_hybrid.fasta")                into ch_sample_hybrid_fasta
    file("${params.name}_refined_high_confidence.tsv") into ch_refined_info_high_conf
  script:
    """
    make_hybrid_database.py \
    --protein_classification $protein_classification \
    --gene_lens $gene_lens \
    --pb_fasta $pb_fasta \
    --gc_fasta $gc_fasta \
    --refined_info $refined_info \
    --pb_cds_gtf $sample_cds \
    --name ${params.name} \
    --lower_kb ${params.lower_kb} \
    --upper_kb ${params.upper_kb} \
    --lower_cpm ${params.lower_cpm} \
    """

}
ch_high_confidence_cds.into{
  ch_high_confidence_cds_track_viz
  ch_high_confidence_cds_multiregion
}

ch_sample_hybrid_fasta.into{
  ch_sample_hybrid_fasta_normal
  ch_sample_hybrid_fasta_rescue
  ch_sample_hybrid_fasta_track_viz
  ch_sample_hybrid_fasta_pep_analysis
}


/*--------------------------------------------------
Mass Spec File Conversion
 * Convert MS .raw files into .mzml files
 * by not specifying the output file, the input file is overwritten
---------------------------------------------------*/
process mass_spec_raw_convert{

    publishDir "${params.outdir}/${params.name}/raw_convert/", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(raw_file) from ch_mass_spec_raw
    output:
        file("*") into ch_mass_spec_converted
    script:
        """
        wine msconvert $raw_file --filter "peakPicking true 1-"
        """
}

ch_uniprot_protein_fasta_uncompressed.into {
  ch_uniprot_protein_fasta_for_metamorphisis
  ch_uniprot_protein_fasta_for_accession_mapping
  ch_uniprot_protein_fasta_for_peptide_novelty
}

ch_mass_spec_combined = ch_mass_spec_mzml.concat(ch_mass_spec_converted)
ch_mass_spec_combined.into{
  ch_mass_spec_for_pacbio_refined
  ch_mass_spec_for_pacbio_filtered
  ch_mass_spec_for_pacbio_hybrid
  ch_mass_spec_for_gencode
  ch_mass_spec_for_uniprot
  ch_mass_spec_for_pacbio_rescue_resolve
}



/*--------------------------------------------------
Metamorpheus GENCODE 
 * Runs Metamorpheus MS search using the GENCODE database
---------------------------------------------------*/
process metamorpheus_with_gencode_database{
    label 'metamorpheus'
    tag "${mass_spec}"
    
    cpus params.max_cpus
    
    publishDir "${params.outdir}/${params.name}/metamorpheus/gencode", mode: 'copy'
    
    when:
      params.mass_spec != false

    input:
        file(gencode_fasta) from ch_gencode_protein_fasta_metamorpheus
        file(mass_spec) from ch_mass_spec_for_gencode.collect()
        file(toml_file) from ch_metamorpheus_toml_gencode
    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv") into ch_gencode_peptides
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv") into ch_gencode_protein_groups
    
    script:
        def toml = toml_file.name != 'NO_TOML_FILE' ? "$toml_file" : 'toml/SearchTask.toml'
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $gencode_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv
        """
}

/*--------------------------------------------------
Metamorpheus UniProt
 * Runs Metamorpheus MS search using the UniProt database
---------------------------------------------------*/
process metamorpheus_with_uniprot_database{
    label 'metamorpheus'
    tag "${mass_spec}"
    cpus params.max_cpus
    publishDir "${params.outdir}/${params.name}/metamorpheus/uniprot", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(toml_file) from ch_metamorpheus_toml_uniprot
        file(uniprot_fasta) from ch_uniprot_protein_fasta_for_metamorphisis
        file(mass_spec) from ch_mass_spec_for_uniprot.collect()

    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv") into ch_uniprot_peptides
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv") into ch_uniprot_protein_groups
    
    script:
        def toml = toml_file.name != 'NO_TOML_FILE' ? "$toml_file" : 'toml/SearchTask.toml'
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $uniprot_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv
        """
}


/*--------------------------------------------------
MetaMorpheus with Sample Specific Database - Refined
 * Runs Metamorpheus MS search using the PacBio refined database
---------------------------------------------------*/
process metamorpheus_with_sample_specific_database_refined{
    label 'metamorpheus'
    tag "${mass_spec}"
    cpus params.max_cpus
    publishDir "${params.outdir}/${params.name}/metamorpheus/pacbio/refined", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(orf_fasta) from ch_refined_fasta_metamorpheus
        file(mass_spec) from ch_mass_spec_for_pacbio_refined.collect()
        file(toml_file) from ch_metamorpheus_toml_pacbio_refined


    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.${params.name}.refined.psmtsv") into ch_pacbio_peptides_refined
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.refined.tsv") into ch_pacbio_protein_groups_refined
    
    script:
        def toml = toml_file.name != 'NO_TOML_FILE' ? "$toml_file" : 'toml/SearchTask.toml'
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $orf_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.refined.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.refined.tsv
        """
}

ch_pacbio_peptides_refined.into{
ch_pacbio_peptides_refined_novel
ch_pacbio_peptides_refined_track_viz
}

/*--------------------------------------------------
MetaMorpheus with Sample Specific Database - Filtered
 * Runs Metamorpheus MS search using the PacBio filtered database
---------------------------------------------------*/
process metamorpheus_with_sample_specific_database_filtered{
    label 'metamorpheus'
    tag "${mass_spec}"
    cpus params.max_cpus
    publishDir "${params.outdir}/${params.name}/metamorpheus/pacbio/filtered", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(orf_fasta) from ch_filtered_protein_fasta_metamorpheus
        file(mass_spec) from ch_mass_spec_for_pacbio_filtered.collect()
        file(toml_file) from ch_metamorpheus_toml_pacbio_filtered


    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.${params.name}.filtered.psmtsv") into ch_pacbio_peptides_filtered
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.filtered.tsv") into ch_pacbio_protein_groups_filtered
    
    script:
        def toml = toml_file.name != 'NO_TOML_FILE' ? "$toml_file" : 'toml/SearchTask.toml'
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $orf_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.filtered.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.filtered.tsv
        """
}

ch_pacbio_peptides_filtered.into{
ch_pacbio_peptides_filtered_novel
ch_pacbio_peptides_filtered_track_viz
}

/*--------------------------------------------------
MetaMorpheus with Sample Specific Database - Hybrid 
 * Runs Metamorpheus MS search using the PacBio hybrid database
---------------------------------------------------*/
process metamorpheus_with_sample_specific_database_hybrid{
    label 'metamorpheus'
    tag "${mass_spec}"
    
    cpus params.max_cpus
    
    publishDir "${params.outdir}/${params.name}/metamorpheus/pacbio/hybrid", mode: 'copy'
    
    when:
      params.mass_spec != false

    input:
        file(orf_fasta) from ch_sample_hybrid_fasta_normal
        file(mass_spec) from ch_mass_spec_for_pacbio_hybrid.collect()
        file(toml_file) from ch_metamorpheus_toml_pacbio_hybrid


    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.${params.name}.hybrid.psmtsv") into ch_pacbio_peptides_hybrid
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.hybrid.tsv") into ch_pacbio_protein_groups_hybrid
    
    script:
        def toml = toml_file.name != 'NO_TOML_FILE' ? "$toml_file" : 'toml/SearchTask.toml'
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $orf_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.hybrid.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.hybrid.tsv
        """
}

ch_pacbio_peptides_hybrid.into{
  ch_pacbio_peptides_hybrid_gtf
  ch_pacbio_peptides_hybrid_novel
  ch_pacbio_peptides_hybrid_track_viz
}

/*--------------------------------------------------
MetaMorpheus using "Rescue and Resolve" Algorithm
 * Runs Metamorpheus MS search using the hybrid database
 * and the version of Metamorpheus containing the
 * "Rescue and Resolve" algorithm. This algorithm allows for
 * "Rescue" of protein isoforms that are normally discarded
 * during protein parsimony, given that the protein isoform
 * has strong evidence of transcriptional expression.
---------------------------------------------------*/
process metamorpheus_with_sample_specific_database_rescue_resolve{
    label 'metamorpheus'
    tag " $mass_spec $orf_fasta $orf_meta  $toml"
    publishDir "${params.outdir}/${params.name}/metamorpheus/pacbio/rescue_resolve", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        // file(orf_calls) from ch_orf_calls
        file(orf_fasta) from ch_sample_hybrid_fasta_rescue
        file(mass_spec) from ch_mass_spec_for_pacbio_rescue_resolve.collect()
        file(toml) from ch_rr_toml
        file(orf_meta) from ch_refined_info_high_conf

    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.${params.name}.rescue_resolve.psmtsv")
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.rescue_resolve.tsv")
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings settings 
        dotnet /metamorpheus/CMD.dll -d $orf_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results --orf $orf_meta --cpm 25
        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.rescue_resolve.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.rescue_resolve.tsv
        """
}

/*--------------------------------------------------
Peptide Analysis
 * Generate a table comparing MS peptide results
 * between the PacBio and GENCODE databases. 
---------------------------------------------------*/
process peptide_analysis{
  publishDir "${params.outdir}/${params.name}/peptide_analysis/", mode: 'copy'
    when:
      params.mass_spec != false
    
    input:
      file(gencode_peptides) from ch_gencode_peptides
      file(gene_isoname) from ch_gene_isoname_pep_analysis
      file(refined_fasta) from ch_refined_fasta_pep_analysis
      file(filtered_fasta) from ch_filtered_protein_fasta_pep_analysis
      file(hybrid_fasta) from ch_sample_hybrid_fasta_pep_analysis
      file(pb_gene) from ch_pb_gene_peptide_analysis

    output:
      file("*")

    script:
      """
      peptide_analysis.py \
      -gmap $gene_isoname \
      --gencode_peptides $gencode_peptides \
      --pb_refined_fasta $refined_fasta \
      --pb_filtered_fasta $filtered_fasta \
      --pb_hybrid_fasta $hybrid_fasta \
      --pb_gene $pb_gene \
      -odir ./
      """
}



/*--------------------------------------------------
Reference Track Visualization
 * Creates tracks to use in UCSC Genome Browser for 
 * Gencode database.
---------------------------------------------------*/
process gencode_track_visualization{
    publishDir "${params.outdir}/${params.name}/track_visualization/reference"

  input:
    file(reference_gtf) from ch_gencode_gtf
  output:
    file("gencode_shaded.bed12")
    file("gencode.filtered.gtf") into ch_gencode_filtered_gtf
  script:
    """
    gencode_filter_protein_coding.py \
    --reference_gtf $reference_gtf

    gtfToGenePred gencode.filtered.gtf gencode.filtered.genePred
    genePredToBed gencode.filtered.genePred gencode.filtered.bed12

    gencode_add_rgb_to_bed.py \
    --gencode_bed gencode.filtered.bed12 \
    --rgb 0,0,140 \
    --version V35
  """
}


/*--------------------------------------------------
Protein Track Visualization
 * Creates tracks to use in UCSC Genome Browser for 
 * refined, filtered, and hybrid PacBio databases.
 * Shades tracks based on abundance (CPMs) and protein classification.
---------------------------------------------------*/
process protein_track_visualization{
  publishDir "${params.outdir}/${params.name}/track_visualization/refined/protein", pattern: "*_refined_*"
  publishDir "${params.outdir}/${params.name}/track_visualization/filtered/protein", pattern: "*_filtered_*"
  publishDir "${params.outdir}/${params.name}/track_visualization/hybrid/protein", pattern: "*_hybrid_*"
  input:
    file(refined_cds) from ch_pb_cds_bed
    file(filtered_cds) from ch_filtered_cds_bed
    file(hybrid_cds) from ch_high_confidence_cds_track_viz
  output:
    file("*_shaded_*")
  script:
    """
    #************************************
    # Convert GTF to Bed
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    gtfToGenePred $refined_cds ${params.name}_refined_cds.genePred
    genePredToBed ${params.name}_refined_cds.genePred ${params.name}_refined_cds.bed12
    #--------------------------
    # Filtered
    #--------------------------
    gtfToGenePred $filtered_cds ${params.name}_filtered_cds.genePred
    genePredToBed ${params.name}_filtered_cds.genePred ${params.name}_filtered_cds.bed12

    #--------------------------
    # Hybrid
    #--------------------------
    gtfToGenePred $hybrid_cds ${params.name}_hybrid_cds.genePred
    genePredToBed ${params.name}_hybrid_cds.genePred ${params.name}_hybrid_cds.bed12


    #************************************
    # Add RGB colors
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    track_add_rgb_colors_to_bed.py \
    --name ${params.name}_refined \
    --bed_file ${params.name}_refined_cds.bed12

    #--------------------------
    # Filtered
    #--------------------------
    track_add_rgb_colors_to_bed.py \
    --name ${params.name}_filtered \
    --bed_file ${params.name}_filtered_cds.bed12

    #--------------------------
    # Hybrid
    #--------------------------
    track_add_rgb_colors_to_bed.py \
    --name ${params.name}_hybrid \
    --bed_file ${params.name}_hybrid_cds.bed12
    """

}






/*--------------------------------------------------
Multiregion BED generation
 * Makes multiregion BED file for UCSC genome browser for 
 * refined, filtered, and hybrid databases. Regions in the
 * BED correpond to coding regions, thereby allowing for
 * intronic regions to be minimized for easier isoform viewing.
---------------------------------------------------*/
process make_multiregion{
  publishDir "${params.outdir}/${params.name}/track_visualization/refined", pattern: "*_refined*"
  publishDir "${params.outdir}/${params.name}/track_visualization/filtered", pattern: "*_filtered*"
  publishDir "${params.outdir}/${params.name}/track_visualization/hybrid", pattern: "*_high_confidence*"
  input:
    file(refined_gtf) from ch_pb_cds_multiregion
    file(filtered_gtf) from ch_filtered_cds_multiregion
    file(reference_gtf) from ch_gencode_filtered_gtf
    file(high_conf_gtf) from ch_high_confidence_cds_multiregion
  output:
    file("*")

  script:
  """
  make_region_bed_for_ucsc.py \
  --name ${params.name}_refined \
  --sample_gtf $refined_gtf \
  --reference_gtf $reference_gtf

  make_region_bed_for_ucsc.py \
  --name ${params.name}_filtered \
  --sample_gtf $filtered_gtf \
  --reference_gtf $reference_gtf

  make_region_bed_for_ucsc.py \
  --name ${params.name}_high_confidence \
  --sample_gtf $high_conf_gtf \
  --reference_gtf $reference_gtf
  """
}

/*--------------------------------------------------
Peptide Track Visualization
 * Makes peptide tracks for UCSC Genome Browser 
 * for refined, filtered and hybrid databases
---------------------------------------------------*/
process peptide_track_visualization{
  publishDir "${params.outdir}/${params.name}/track_visualization/refined/peptide", pattern: "*_refined_*"
  publishDir "${params.outdir}/${params.name}/track_visualization/filtered/peptide", pattern: "*_filtered_*"
  publishDir "${params.outdir}/${params.name}/track_visualization/hybrid/peptide", pattern: "*_hybrid_*"

  when:
    params.mass_spec != false
  
  input:
    file(sample_gtf) from ch_pb_cds_peptide_gtf
    file(reference_gtf) from ch_gencode_gtf
    file(refined_peptides) from ch_pacbio_peptides_refined_track_viz
    file(filtered_peptides) from ch_pacbio_peptides_filtered_track_viz
    file(hybrid_peptides) from ch_pacbio_peptides_hybrid_track_viz
    file(pb_gene) from ch_pb_gene_peptide_gtf
    file(gene_isoname) from ch_gene_isoname_pep_viz
    file(hybrid_fasta) from ch_sample_hybrid_fasta_track_viz
    file(refined_fasta) from ch_refined_fasta_peptide_viz
  output:
    file("*.gtf")
    file("*.bed12")
  script:
    """
    #************************************
    # Make peptide gtf files
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    make_peptide_gtf_file.py \
    --name ${params.name}_refined \
    --sample_gtf $sample_gtf \
    --reference_gtf $reference_gtf \
    --peptides $refined_peptides \
    --pb_gene $pb_gene \
    --gene_isoname $gene_isoname \
    --refined_fasta $refined_fasta

    #--------------------------
    # Filtered
    #--------------------------
    make_peptide_gtf_file.py \
    --name ${params.name}_filtered \
    --sample_gtf $sample_gtf \
    --reference_gtf $reference_gtf \
    --peptides $filtered_peptides \
    --pb_gene $pb_gene \
    --gene_isoname $gene_isoname \
    --refined_fasta $refined_fasta

    #--------------------------
    # Hybrid
    #--------------------------
    make_peptide_gtf_file.py \
    --name ${params.name}_hybrid \
    --sample_gtf $sample_gtf \
    --reference_gtf $reference_gtf \
    --peptides $hybrid_peptides \
    --pb_gene $pb_gene \
    --gene_isoname $gene_isoname \
    --refined_fasta $hybrid_fasta

    #************************************
    # Convert GTF to bed and add RGB
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    gtfToGenePred ${params.name}_refined_peptides.gtf ${params.name}_refined_peptides.genePred
    genePredToBed ${params.name}_refined_peptides.genePred ${params.name}_refined_peptides.bed12
    # add rgb to colorize specific peptides 
    finalize_peptide_bed.py \
    --bed ${params.name}_refined_peptides.bed12 \
    --name ${params.name}_refined

    #--------------------------
    # Filtered
    #--------------------------
    gtfToGenePred ${params.name}_filtered_peptides.gtf ${params.name}_filtered_peptides.genePred
    genePredToBed ${params.name}_filtered_peptides.genePred ${params.name}_filtered_peptides.bed12
    # add rgb to colorize specific peptides 
    finalize_peptide_bed.py \
    --bed ${params.name}_filtered_peptides.bed12 \
    --name ${params.name}_filtered
    #--------------------------
    # High Confidence
    #--------------------------
    gtfToGenePred ${params.name}_hybrid_peptides.gtf ${params.name}_hybrid_peptides.genePred
    genePredToBed ${params.name}_hybrid_peptides.genePred ${params.name}_hybrid_peptides.bed12
    # add rgb to colorize specific peptides 
    finalize_peptide_bed.py \
    --bed ${params.name}_hybrid_peptides.bed12 \
    --name ${params.name}_hybrid
    """
  
}



/*--------------------------------------------------
Accession Mapping 
 * Maps protein entries within GENCODE, UniProt and PacBio databases to one
 * another based on sequence similarity 
---------------------------------------------------*/
process accession_mapping{
  publishDir "${params.outdir}/${params.name}/accession_mapping/", mode: 'copy'
  when:
    params.mass_spec != false

  input:
    file(pacbio_fasta) from ch_refined_fasta_mapping
    file(gencode_fasta) from ch_gencode_protein_fasta_mapping
    file(uniprot_fasta) from ch_uniprot_protein_fasta_for_accession_mapping
  
  output:
    file("accession_map_gencode_uniprot_pacbio.tsv") into ch_accession_map
    file("*")
  
  script:
    """
    accession_mapping.py \
    --gencode_fasta $gencode_fasta \
    --pacbio_fasta $pacbio_fasta \
    --uniprot_fasta $uniprot_fasta \
    """
}


/*--------------------------------------------------
Protein Group Comparison
 * Determine the relationship between protein groups identified in Metamorpheus 
 * using GENCODE, UniProt, and/or PacBio databases
---------------------------------------------------*/
process protein_group_compare{
    publishDir "${params.outdir}/${params.name}/protein_group_compare/", mode: 'copy'
    when:
      params.mass_spec != false
    input: 
      file(pacbio_protein_groups) from ch_pacbio_protein_groups_hybrid
      file(gencode_protein_groups) from ch_gencode_protein_groups
      file(uniprot_protein_groups) from ch_uniprot_protein_groups
      file(mapping) from ch_accession_map
    output:
      file("*")
    script:
      """
      protein_groups_compare.py \
      --pg_fileOne $gencode_protein_groups \
      --pg_fileTwo $pacbio_protein_groups \
      --mapping $mapping \
      --output ./

      protein_groups_compare.py \
      --pg_fileOne $gencode_protein_groups \
      --pg_fileTwo $uniprot_protein_groups \
      --mapping $mapping \
      --output ./

      protein_groups_compare.py \
      --pg_fileOne $uniprot_protein_groups \
      --pg_fileTwo $pacbio_protein_groups \
      --mapping $mapping \
      --output ./
      """
}

/*--------------------------------------------------
Novel Peptides
 * Finds novel peptides between sample database and GENCODE 
 * Novel peptide is defined as a peptide found in PacBio that 
 * could not be found in Gencode
 * Novel peptides found for refined, filtered, and hybrid databases
---------------------------------------------------*/
process peptide_novelty_analysis{
  publishDir "${params.outdir}/${params.name}/novel_peptides/", mode: 'copy'
  input:
    file(peptides_refined) from ch_pacbio_peptides_refined_novel
    file(peptides_filtered) from ch_pacbio_peptides_filtered_novel
    file(peptides_hybrid) from ch_pacbio_peptides_hybrid_novel
    file(gencode_fasta) from ch_gencode_protein_fasta_novel
    file(uniprot_fasta) from ch_uniprot_protein_fasta_for_peptide_novelty
  output:
    file("*")
  
  script:
    """
    #--------------------------
    # Refined
    #--------------------------
    peptide_novelty_analysis.py \
    --pacbio_peptides $peptides_refined \
    --gencode_fasta $gencode_fasta \
    --uniprot_fasta $uniprot_fasta \
    --name ${params.name}_refined

    #--------------------------
    # Filtered
    #--------------------------
    peptide_novelty_analysis.py \
    --pacbio_peptides $peptides_filtered \
    --gencode_fasta $gencode_fasta \
    --uniprot_fasta $uniprot_fasta \
    --name ${params.name}_filtered

    #--------------------------
    # High Confidence
    #--------------------------
    peptide_novelty_analysis.py \
    --pacbio_peptides $peptides_hybrid \
    --gencode_fasta $gencode_fasta \
    --uniprot_fasta $uniprot_fasta \
    --name ${params.name}_hybrid
    """

}


def logHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
    ${c_cyan}  sheynkman-lab/Long-Read-Proteogenomics v${workflow.manifest.version}${c_reset}
    ${c_cyan}███████╗██╗  ██╗███████╗██╗   ██╗███╗   ██╗██╗  ██╗███╗   ███╗ █████╗ ███╗   ██╗      ██╗      █████╗ ██████╗            ${c_reset}       
    ${c_cyan}██╔════╝██║  ██║██╔════╝╚██╗ ██╔╝████╗  ██║██║ ██╔╝████╗ ████║██╔══██╗████╗  ██║      ██║     ██╔══██╗██╔══██╗           ${c_reset}        
    ${c_cyan}███████╗███████║█████╗   ╚████╔╝ ██╔██╗ ██║█████╔╝ ██╔████╔██║███████║██╔██╗ ██║█████╗██║     ███████║██████╔╝           ${c_reset}
    ${c_cyan}╚════██║██╔══██║██╔══╝    ╚██╔╝  ██║╚██╗██║██╔═██╗ ██║╚██╔╝██║██╔══██║██║╚██╗██║╚════╝██║     ██╔══██║██╔══██╗           ${c_reset}
    ${c_cyan}███████║██║  ██║███████╗   ██║   ██║ ╚████║██║  ██╗██║ ╚═╝ ██║██║  ██║██║ ╚████║      ███████╗██║  ██║██████╔╝           ${c_reset}
    ${c_cyan}╚══════╝╚═╝  ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝      ╚══════╝╚═╝  ╚═╝╚═════╝            ${c_reset}
    ${c_cyan}                                                                                                                         ${c_reset}
    ${c_cyan}██╗      ██████╗ ███╗   ██╗ ██████╗     ██████╗ ███████╗ █████╗ ██████╗                                                  ${c_reset}
    ${c_cyan}██║     ██╔═══██╗████╗  ██║██╔════╝     ██╔══██╗██╔════╝██╔══██╗██╔══██╗                                                 ${c_reset}
    ${c_cyan}██║     ██║   ██║██╔██╗ ██║██║  ███╗    ██████╔╝█████╗  ███████║██║  ██║                                                 ${c_reset}
    ${c_cyan}██║     ██║   ██║██║╚██╗██║██║   ██║    ██╔══██╗██╔══╝  ██╔══██║██║  ██║                                                 ${c_reset}
    ${c_cyan}███████╗╚██████╔╝██║ ╚████║╚██████╔╝    ██║  ██║███████╗██║  ██║██████╔╝                                                 ${c_reset}
    ${c_cyan}╚══════╝ ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═════╝                                                  ${c_reset}
    ${c_cyan}                                                                                                                         ${c_reset}
    ${c_cyan}██████╗ ██████╗  ██████╗ ████████╗███████╗ ██████╗  ██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗██╗ ██████╗███████╗    ${c_reset}
    ${c_cyan}██╔══██╗██╔══██╗██╔═══██╗╚══██╔══╝██╔════╝██╔═══██╗██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██║██╔════╝██╔════╝    ${c_reset}
    ${c_cyan}██████╔╝██████╔╝██║   ██║   ██║   █████╗  ██║   ██║██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║██║██║     ███████╗    ${c_reset}
    ${c_cyan}██╔═══╝ ██╔══██╗██║   ██║   ██║   ██╔══╝  ██║   ██║██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██║██║     ╚════██║    ${c_reset}
    ${c_cyan}██║     ██║  ██║╚██████╔╝   ██║   ███████╗╚██████╔╝╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║██║╚██████╗███████║    ${c_reset}
    ${c_cyan}╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚══════╝ ╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝    ${c_reset}
    ${c_cyan}                                                                                                                                       
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}





