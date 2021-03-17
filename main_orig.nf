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
// Header log info
log.info "\nPARAMETERS SUMMARY"
log.info "mainScript                            : ${params.mainScript}"
log.info "config                                : ${params.config}"
log.info "max_cpus                              : ${params.max_cpus}"
log.info "outdir                                : ${params.outdir}"
log.info "name                                  : ${params.name}"
log.info "gencode_gtf                           : ${params.gencode_gtf}"
log.info "gencode_transcript_fasta              : ${params.gencode_transcript_fasta}"
log.info "gencode_translation_fasta             : ${params.gencode_translation_fasta}"
log.info "genome_fasta                          : ${params.genome_fasta}"
log.info "fastq_read_1                          : ${params.fastq_read_1}"
log.info "fastq_read_2                          : ${params.fastq_read_2}"
log.info "star_genome_dir                       : ${params.star_genome_dir}"
log.info "sample_ccs                            : ${params.sample_ccs}"
log.info "primers_fasta                         : ${params.primers_fasta}"
log.info "hexamer                               : ${params.hexamer}"
log.info "logit_model                           : ${params.logit_model}"
log.info "sample_kallisto_tpm                   : ${params.sample_kallisto_tpm}"
log.info "normalized_ribo_kallisto              : ${params.normalized_ribo_kallisto}"
log.info "uniprot_fasta                         : ${params.uniprot_fasta}"
log.info "uniprot_protein_fasta                 : ${params.uniprot_protein_fasta}"
log.info "protein_coding_only                   : ${params.protein_coding_only}"
log.info "refine_cutoff                         : ${params.refine_cutoff}"
log.info "mass_spec                             : ${params.mass_spec}"
log.info ""

// if (!params.gencode_gtf) exit 1, "Cannot find gtf file for parameter --gencode_gtf: ${params.gencode_gtf}"
ch_gencode_gtf = Channel.value(file(params.gencode_gtf))

if (!params.gencode_transcript_fasta) exit 1, "Cannot find any file for parameter --gencode_transcript_fasta: ${params.gencode_transcript_fasta}"
ch_gencode_transcript_fasta= Channel.value(file(params.gencode_transcript_fasta))

if (!params.gencode_translation_fasta) exit 1, "Cannot find any file for parameter --gencode_translation_fasta: ${params.gencode_translation_fasta}"

if (params.gencode_translation_fasta.endsWith('.gz')){
ch_gencode_translation_fasta = Channel.value(file(params.gencode_translation_fasta))
}
if (!params.gencode_translation_fasta.endsWith('.gz')){
ch_gencode_translation_fasta_uncompressed = Channel.value(file(params.gencode_translation_fasta))
}

if (!params.sample_ccs) exit 1, "Cannot find file for parameter --sample_ccs: ${params.sample_ccs}"
ch_sample_ccs = Channel.value(file(params.sample_ccs))

if (!params.genome_fasta) exit 1, "Cannot find any seq file for parameter --genome_fasta: ${params.genome_fasta}"
ch_genome_fasta = Channel.value(file(params.genome_fasta))

if (!params.primers_fasta) exit 1, "Cannot find any seq file for parameter --primers_fasta: ${params.primers_fasta}"
ch_primers_fasta = Channel.value(file(params.primers_fasta))

if (!params.hexamer) exit 1, "Cannot find headmer file for parameter --hexamer: ${params.hexamer}"
ch_hexamer = Channel.value(file(params.hexamer))

if (!params.logit_model) exit 1, "Cannot find any logit model file for parameter --logit_model: ${params.logit_model}"
ch_logit_model =  Channel.value(file(params.logit_model))

if (!params.sample_kallisto_tpm) exit 1, "Cannot find any sample_kallisto_tpm file for parameter --sample_kallisto_tpm: ${params.sample_kallisto_tpm}"
ch_sample_kallisto = Channel.value(file(params.sample_kallisto_tpm))

if (!params.normalized_ribo_kallisto) exit 1, "Cannot find any normalized_ribo_kallisto file for parameter --normalized_ribo_kallisto: ${params.normalized_ribo_kallisto}"
ch_normalized_ribo_kallisto = Channel.value(file(params.normalized_ribo_kallisto))

if (!params.uniprot_protein_fasta) exit 1, "Cannot find any file for parameter --uniprot_protein_fasta: ${params.uniprot_protein_fasta}"
ch_uniprot_protein_fasta = Channel.value(file(params.uniprot_protein_fasta))

if (!params.fastq_read_1) exit 1, "No file found for the parameter --fastq_read_1 at the location ${params.fastq_read_1}"
if (!params.fastq_read_2) exit 1, "No file found for the parameter --fastq_read_2 at the location ${params.fastq_read_2}"
ch_fastq_reads = Channel.from(params.fastq_read_1, params.fastq_read_2).filter(String).flatMap{ files(it) }


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

ch_genome_fasta.into{
  ch_genome_fasta_star
  ch_genome_fasta_isoseq
  ch_genome_fasta_sqanti
}


/*--------------------------------------------------
Gencode Database
---------------------------------------------------*/

if (params.gencode_translation_fasta.endsWith('.gz')) {
  process gunzip_gencode_translation_fasta {
  tag "decompress gzipped fasta"
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

process make_gencode_database {
  tag "${gencode_translation_fasta}"
  cpus 1
  publishDir "${params.outdir}/gencode_db/", mode: 'copy'

  input:
  file(gencode_translation_fasta) from ch_gencode_translation_fasta_uncompressed
  
  output:
  file("gencode_protein.fasta") into ch_gencode_protein_fasta
  file("gencode_isoname_clusters.tsv")
  
  script:
  """
  make_gencode_database.py \
  --gencode_fasta $gencode_translation_fasta \
  --output_fasta gencode_protein.fasta \
  --output_cluster gencode_isoname_clusters.tsv
  """
}

ch_gencode_protein_fasta.into{
  ch_gencode_protein_fasta_metamorpheus
  ch_gencode_protein_fasta_mapping
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

if(params.star_genome_dir != false){
    Channel
    .fromPath(params.star_genome_dir, type:'dir')
    .set{ch_genome_dir}
}
else{
    process star_generate_genome{
        cpus params.max_cpus
        publishDir "${params.outdir}/star_index", mode: "copy"
        
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
        publishDir "${params.outdir}/star", mode: "copy"
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
---------------------------------------------------*/

process sqanti3 {
  tag "${fl_count}, ${gencode_gtf}, ${gencode_fasta}, ${sample_gtf},"
  cpus params.max_cpus
  publishDir "${params.outdir}/sqanti3/", mode: 'copy'

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
    --minimum_illumina_coverage 3 \
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
---------------------------------------------------*/

process six_frame_translation {
    cpus 1
    tag "${classification}, ${ensg_gene}"
    publishDir "${params.outdir}/pacbio_6frm_gene_grouped/", mode: 'copy'

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
---------------------------------------------------*/

process transcriptome_summary {
  cpus 1
  publishDir "${params.outdir}/transcriptome_summary/", mode: 'copy'

  input:
  file(sqanti_classification) from ch_sample_classification_transcriptome
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

ch_pb_gene.into{
  ch_pb_gene_orf
  ch_pb_gene_cds
  ch_pb_gene_peptide_analysis
  ch_pb_gene_peptide_gtf
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
}

/*--------------------------------------------------
Refined DB Generation 
---------------------------------------------------*/

process refine_orf_database {
  cpus 1
  tag "${best_orfs}, ${sample_fasta}, ${params.protein_coding_only}, ${protein_coding_genes}, ${params.refine_cutoff}" 

  publishDir "${params.outdir}/refined_database/", mode: 'copy'

  input:
  file(best_orfs) from ch_best_orf_refine
  file(sample_fasta) from ch_sample_fasta_refine
  file(protein_coding_genes) from ch_protein_coding_genes_db
  
  output:
  file("*")
  file("${params.name}_orf_refined.tsv") into ch_refined_info
  file("${params.name}_orf_refined.fasta") into ch_refined_fasta
  
  script:
  """
  refine_orf_database.py \
  --name ${params.name} \
  --orfs $best_orfs \
  --pb_fasta $sample_fasta \
  --coding_score_cutoff ${params.refine_cutoff} \
  """
}

ch_refined_fasta.into{
  ch_refined_fasta_metamorpheus
  ch_refined_fasta_pep_analysis
  ch_refined_fasta_peptide_gtf
  ch_refined_fasta_mapping
}




if(params.mass_spec != false){
  Channel
    .fromPath("${params.mass_spec}/*.raw")
    .set{ch_mass_spec_raw}
  Channel
      .fromPath("${params.mass_spec}/*.{mzml,mzML}")
      .set{ch_mass_spec_mzml}
}
else{
  Channel
    .from("no mass spec")
    .set{ch_mass_spec_raw}
  Channel
    .from("no mass spec")
    .set{ch_mass_spec_mzml}

}


/*--------------------------------------------------
MetaMorpheus wtih Sample Specific Database
---------------------------------------------------*/

process mass_spec_raw_convert{
    // publishDir "${params.outdir}/raw_convert/", mode: 'copy'
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
ch_mass_spec_combined.into{
  ch_mass_spec_for_pacbio
  ch_mass_spec_for_gencode
  ch_mass_spec_for_uniprot
  ch_mass_spec_for_pacbio_rescue_resolve
}

process metamorpheus_with_sample_specific_database{
    tag "${mass_spec}"
    cpus params.max_cpus
    publishDir "${params.outdir}/metamorpheus/pacbio", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(orf_fasta) from ch_refined_fasta_metamorpheus
        file(mass_spec) from ch_mass_spec_for_pacbio.collect()

    output:
        file("toml/*")
        file("search_results/Task1SearchTask/All*")
        file("search_results/Task1SearchTask/prose.txt")
        file("search_results/Task1SearchTask/results.txt")
        file("search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv") into ch_pacbio_peptides
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv") into ch_pacbio_protein_groups
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $orf_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.tsv
        """
}

process metamorpheus_with_gencode_database{
    tag "${mass_spec}"
    cpus params.max_cpus
    publishDir "${params.outdir}/metamorpheus/gencode", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(gencode_fasta) from ch_gencode_protein_fasta_metamorpheus
        file(mass_spec) from ch_mass_spec_for_gencode.collect()

    output:
        file("toml/*")
        file("search_results/Task1SearchTask/All*")
        file("search_results/Task1SearchTask/prose.txt")
        file("search_results/Task1SearchTask/results.txt")
        file("search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv") into ch_gencode_peptides
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv") into ch_gencode_protein_groups
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $gencode_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.tsv
        """
}

process metamorpheus_with_uniprot_database{
    tag "${mass_spec}"
    cpus params.max_cpus
    publishDir "${params.outdir}/metamorpheus/uniprot", mode: 'copy'
    when:
      params.mass_spec != false

    input:
        file(uniprot_fasta) from ch_uniprot_protein_fasta
        file(mass_spec) from ch_mass_spec_for_uniprot.collect()

    output:
        file("toml/*")
        file("search_results/Task1SearchTask/All*")
        file("search_results/Task1SearchTask/prose.txt")
        file("search_results/Task1SearchTask/results.txt")
        file("search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv") into ch_uniprot_peptides
        file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv") into ch_uniprot_protein_groups
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
        dotnet /metamorpheus/CMD.dll -d $uniprot_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results

        mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv
        mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv
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
      file(refined_fasta) from ch_refined_fasta_pep_analysis
      file(six_frame) from ch_six_frame
      file(pb_gene) from ch_pb_gene_peptide_analysis
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
PacBio CDS GTF 
---------------------------------------------------*/

process make_pacbio_cds_gtf {
  cpus 1

  publishDir "${params.outdir}/pacbio_cds/", mode: 'copy'

  input:
    file(sample_gtf) from ch_sample_gtf_cds
    file(refined_info) from ch_refined_info
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


/*--------------------------------------------------
Novel Peptides
---------------------------------------------------*/
// TODO - proposed module to list all peptides from sample search and mark novel peptides


/*--------------------------------------------------
Accession Mapping 
---------------------------------------------------*/

process accession_mapping{
  publishDir "${params.outdir}/accession_mapping/", mode: 'copy'
  when:
    params.mass_spec != false

  input:
    file(pacbio_fasta) from ch_refined_fasta_mapping
    file(gencode_fasta) from ch_gencode_protein_fasta_mapping
    file(uniprot_fasta) from ch_uniprot_protein_fasta
  
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
---------------------------------------------------*/
process protein_group_compare{
    publishDir "${params.outdir}/protein_group_compare/", mode: 'copy'
    when:
      params.mass_spec != false
    input: 
      file(pacbio_protein_groups) from ch_pacbio_protein_groups
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






