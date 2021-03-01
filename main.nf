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
 * Anne Deslattes Mays (adeslat@scitechcon.org)
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
    .ifEmpty { error "Cannot find any gencode_fasta file for parameter --gencode_transcript_fasta: ${params.gencode_transcript_fasta}" }
    .set { ch_gencode_transcript_fasta }  

Channel
    .value(file(params.gencode_translation_fasta))
    .ifEmpty { error "Cannot find any gencode_fasta file for parameter --gencode_translation_fasta: ${params.gencode_translation_fasta}" }
    .set { ch_gencode_translation_fasta }  

Channel
    .value(file(params.sample_ccs))
    .ifEmpty { error "Cannot find file for parameter --sample_ccs: ${params.sample_ccs}" }
    .set { ch_sample_ccs }   

// TODO - rename this to "--genome_fasta"
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

Channel
    .from(params.fastq_read_1, params.fastq_read_2)
    .filter(String)
    .flatMap{ files(it) }
    .set{ ch_fastq_reads }


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

// partition channels for use by multiple modules
ch_protein_coding_genes.into{
  ch_protein_coding_genes_db
  ch_protein_coding_genes_filter
}

ch_ensg_gene.into{
  ch_ensg_gene_filter
  ch_ensg_gene_six_frame
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
  // TODO - change to ch_gencode_protein_fasta
  file("gencode_protein.fasta") into ch_genome_protein_fasta
  // TODO - what happens when a file doesn't go into a variable - does it get output to results?
  file("gencode_isoname_clusters.tsv")
  
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
  // TODO - rename ch_genome_fasta
  file(gencode_fasta) from ch_gencode_fasta
  file(primers_fasta) from ch_primers_fasta
  
  output:
  // TODO - decide which files to keep in results folder, asked Liz about this
  file("*")
  file("${params.name}.collapsed.gff") into ch_isoseq_gtf
  // TODO - no longer needed because getting fasta from sqanti?
  // file("${params.name}.collapsed.fasta") into ch_isoseq_fasta
  file("${params.name}.collapsed.abundance.txt") into ch_fl_count
  file("${params.name}")
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
  pbmm2 align $gencode_fasta ${params.name}.clustered.hq.bam ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${task.cpus} --log-level INFO

  # collapse redundant reads
  isoseq3 collapse ${params.name}.aligned.bam ${params.name}.collapsed.gff
  """
}




// TODO - what was this code snippet for? runnning partially?
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
STAR Alignment
---------------------------------------------------*/

// TODO - "== true" ?   
if(params.star_genome_dir != false){
    Channel
    .fromPath(params.star_genome_dir, type:'dir')
    .set{ch_genome_dir}
}
else{
    process star_generate_genome{
        cpus params.max_cpus
        
        when:
        (params.fastq_read_1 != false | params.fastq_read_2 !=false) & params.star_genome_dir == false

        input :
            file(gencode_gtf) from ch_gencode_gtf
            file(gencode_fasta) from ch_gencode_fasta

        output:
            path("star_genome") into ch_genome_dir

        script:
        """
        mkdir star_genome
        STAR --runThreadN  ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir star_genome \
        --genomeFastaFiles $gencode_fasta \
        --sjdbGTFfile $gencode_gtf \
        --genomeSAindexNbases 11
        """
    }
}

if(params.fastq_read_1 != false | params.fastq_read_2 !=false){
    process star_alignment{
        cpus params.max_cpus
        publishDir "${params.outdir}/star", mode: "copy"
        when:
            params.fastq_read_1 != false | params.fastq_read_2 !=false

        input :
            file(fastq_reads) from ch_fastq_reads.collect()
            path(genome_dir) from ch_genome_dir

        output:
            // TODO - don't recommend keeping all files, the BAM/SAM files are large
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
---------------------------------------------------*/

process sqanti3 {
  tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"
  cpus params.max_cpus
  publishDir "${params.outdir}/sqanti3/", mode: 'copy'

  input:
  file(fl_count) from ch_fl_count
  file(gencode_gtf) from ch_gencode_gtf
  // TODO - change to ch_genome_fasta, also need to change in commands below
  file(gencode_fasta) from ch_gencode_fasta
  file(sample_gtf) from ch_isoseq_gtf
  file(star_junction) from ch_star_junction
  
  
  output:
  file("${params.name}_classification.txt") into ch_sample_unfiltered_classification
  file("${params.name}_corrected.fasta") into ch_sample_unfiltered_fasta
  file("${params.name}_corrected.gtf") into ch_sample_unfiltered_gtf
  file("*")
  
  script:
  if(params.fastq_read_1 == false & params.fastq_read_2 ==false)
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
  else
    """
    sqanti3_qc.py \
    $sample_gtf \
    $gencode_gtf \
    $gencode_fasta \
    --skipORF \
    -o ${params.name} \
    --fl_count $fl_count  \
    --gtf \
    -c $star_junction
    """
  //
}

// TODO - add an additional filter for only nnc with illumina coverage
process filter_sqanti {
  publishDir "${params.outdir}/sqanti3-filtered/", mode: 'copy'

  input:
    file(classification) from ch_sample_unfiltered_classification
    file(sample_fasta) from ch_sample_unfiltered_fasta
    file(sample_gtf) from ch_sample_unfiltered_gtf
    file(protein_coding_genes) from ch_protein_coding_genes
    file(ensg_gene) from ch_ensg_gene_filter

  output:
    file("filtered_${params.name}_classification.txt") into ch_sample_classification
    file("filtered_${params.name}_corrected.fasta") into ch_sample_fasta
    file("filtered_${params.name}_corrected.gtf") into ch_sample_gtf

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
    """
}

// TODO - does this need to stay?
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
    file(ensg_gene) from ch_ensg_gene_six_frame
    file(sample_fasta) from ch_sample_fasta

    output:
    file("${params.name}.6frame.fasta") into ch_six_frame

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
  
  // TODO - sqanti_isoform_info.tsv outputs isonames with underscores, convert back to hyphens
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

// TODO - why put into different channel names versus using the same channel?
ch_pb_gene.into{
  ch_pb_gene_orf
  ch_pb_gene_cds
  ch_pb_gene_peptide
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

  // TODO - do we use "jurkat.ORF_seqs.fa"? downstream, potentially not output? 
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


/*--------------------------------------------------
ORF Calling 
---------------------------------------------------*/

process orf_calling {
  tag "${orf_coord}, ${gencode_gtf}, ${sample_gtf}, ${pb_gene}, ${classification}, ${sample_fasta} "
  cpus params.max_cpus
  publishDir "${params.outdir}/orf_calling/", mode: 'copy'

  input:
  file(cpat_orfs) from ch_cpat_all_orfs_for_orf_calling
  file(gencode_gtf) from ch_gencode_gtf
  file(sample_gtf) from ch_sample_gtf
  file(sample_fasta) from ch_sample_fasta
  file(pb_gene) from ch_pb_gene_orf
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
  // TODO - would short params go here too? or do they get called within the script line (as it is now)
  file(best_orfs) from ch_best_orf
  file(sample_fasta) from ch_sample_fasta
  file(protein_coding_genes) from ch_protein_coding_genes_db
  
  output:
  // TODO - Ben will modify python script so output is only <sample>_orf_refined.tsv (previously *_orf_aggregated.tsv), same for fasta
  file("*")
  file("${params.name}_orf_aggregated.tsv") into ch_refined_info
  file("${params.name}_orf_aggregated.fasta") into ch_refined_fasta
  
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
// TODO - complete the module?
process make_pacbio_cds_gtf {
  cpus 1

  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
  file(sample_gtf) from ch_sample_gtf
  file(agg_orfs) from ch_refined_info
  file(refined_orfs) from ch_best_orf
  file(pb_gene) from ch_pb_gene_cds
  
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


if params.mass_spec != false{
  Channel
    .fromPath("${params.mass_spec}/*.raw")
    .set{ch_mass_spec_raw}

  Channel
      .fromPath("${params.mass_spec}/*.{mzml,mzML}")
      .set{ch_mass_spec_mzml}
}


/*--------------------------------------------------
MetaMorpheus wtih Sample Specific Database
---------------------------------------------------*/

process mass_spec_raw_convert{
    publishDir "${params.outdir}/raw_convert/", mode: 'copy'
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

ch_mass_spec_combined = ch_mass_spec_mzml.concat(ch_mass_spec_converted)

process metamorpheus_with_sample_specific_database{
    tag "${mass_spec}"
    publishDir "${params.outdir}/metamorpheus/", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(orf_fasta) from ch_orf_fasta
        file(mass_spec) from ch_mass_spec_combined.collect()

    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.psmtsv") into ch_sample_specific_peptides
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings settings
        dotnet /metamorpheus/CMD.dll -d $orf_fasta -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
        """
}

process metamorpheus_with_gencode_database{
    tag "${mass_spec}"
    publishDir "${params.outdir}/metamorpheus/", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(gencode_fasta) from ch_genome_protein_fasta
        file(mass_spec) from ch_mass_spec_combined.collect()

    output:
        file("toml/*")
        file("search_results/*")
        file("search_results/Task1SearchTask/AllPeptides.psmtsv") into ch_gencode_peptides
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings settings
        dotnet /metamorpheus/CMD.dll -d $gencode_fasta -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
        """
}

/*--------------------------------------------------
Peptide Analysis
---------------------------------------------------*/

process peptide_analysis{
  publishDir "${params.outdir}/peptide_analysis/", mode: 'copy'
    when:
      params.mass_spec != false
    
    input:
      file(gencode_peptides) from ch_gencode_peptides
      file(gene_isoname) from ch_gene_isoname
      file(refined_fasta) from ch_refined_fasta
      file(six_frame) from ch_six_frame
      file(pb_gene) from ch_pb_gene_peptide
      file(cpat_all_orfs) from ch_cpat_all_orfs_for_peptide_analysis
      file(cpat_best_orf) from ch_cpat_best_orf
      file(cpat_protein_fasta) from ch_cpat_protein_fasta

    output:
      file("*")

    script:
      """
      peptide_analysis.py \
      -gmap $gene_isoname \
      -gc $gencode_peptides \
      -pb $refined_fasta \
      -sft $six_frame \
      --pb_gene $pb_gene \
      --cpat_all_orfs $cpat_all_orfs \
      --cpat_best_orf $cpat_best_orf \
      --cpat_orf_protein_fasta $cpat_protein_fasta \
      -odir ./
      """
}



/*--------------------------------------------------
Novel Peptides
---------------------------------------------------*/
// TODO - proposed module to list all peptides from sample search and mark novel peptides


/*--------------------------------------------------
Accession Mapping 
---------------------------------------------------*/
// Have "exact" sequence accession mapping (simple py script)
// Other option is blast-based or fuzzy-matching-based sequence comparisons


/*--------------------------------------------------
Visualization Tracks
---------------------------------------------------*/
// TODO - implement this module
// make cds gtf (already done in other module)
// make peptide gtf
// make custom ranges -> option to do this for:
//  1) gencode + pacbio cds ranges
//  2) gencode + pacbio exon/cds ranges
//  3) pacbio cds ranges only
//  4) pacbio exon/cds ranges
// additional options for:
//  1) color by abundance
//  2) color peptides by uniqueness
//  3) color pacbio transcripts by novelty


/*--------------------------------------------------
Protein Inference Analysis
---------------------------------------------------*/
// TODO - implement Rachel's code that does a cross-comparison of protein groups
// NOTE - her code requires a map from the accession mapping module


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






