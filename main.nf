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
log.info "sample_ccs                            : ${params.sample_ccs}"
log.info "primers_fasta                         : ${params.primers_fasta}"
log.info ""


if (params.gencode_gtf.endsWith('.gz')){
ch_gencode_gtf = Channel.value(file(params.gencode_gtf))
}
if (!params.gencode_gtf.endsWith('.gz')){
ch_gencode_gtf_uncompressed = Channel.value(file(params.gencode_gtf))
}

if (params.genome_fasta.endsWith('.gz')){
ch_genome_fasta = Channel.value(file(params.genome_fasta))
}
if (!params.genome_fasta.endsWith('.gz')){
ch_genome_fasta_uncompressed = Channel.value(file(params.genome_fasta))
}

if (!params.sample_ccs) exit 1, "Cannot find file for parameter --sample_ccs: ${params.sample_ccs}"
ch_sample_ccs = Channel.value(file(params.sample_ccs))

if (!params.primers_fasta) exit 1, "Cannot find any seq file for parameter --primers_fasta: ${params.primers_fasta}"
ch_primers_fasta = Channel.value(file(params.primers_fasta))


if (params.genome_fasta.endsWith('.gz')) {
  process gunzip_gencode_fasta {
  tag "decompress gzipped fasta"
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



if (params.gencode_gtf.endsWith('.gz')) {
  process gunzip_gencode_gtf {
  tag "decompress gzipped gtf"
  cpus 1

  input:
  file(gencode_gtf) from ch_gencode_gtf

  output:
  file("*.gtf") into ch_gencode_gtf_uncompressed

  script:
  """
  gunzip -f ${gencode_gtf}
  """
  }
}

ch_genome_fasta_uncompressed.into{
    ch_genome_fasta_isoseq
    ch_genome_fasta_sqanti
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


/*--------------------------------------------------
SQANTI3
---------------------------------------------------*/

process sqanti3 {
  tag "${fl_count}, ${gencode_gtf}, ${genome_fasta}, ${sample_gtf},"
  cpus params.max_cpus
  publishDir "${params.outdir}/sqanti3/", mode: 'copy'

  input:
  file(fl_count) from ch_fl_count
  file(gencode_gtf) from ch_gencode_gtf_uncompressed
  file(genome_fasta) from ch_genome_fasta_sqanti
  file(sample_gtf) from ch_isoseq_gtf
  
  
  output:
  file("${params.name}_classification.txt") into ch_sample_unfiltered_classification
  file("${params.name}_corrected.fasta") into ch_sample_unfiltered_fasta
  file("${params.name}_corrected.gtf") into ch_sample_unfiltered_gtf
  file("${params.name}_junctions.txt")
  file("${params.name}_sqanti_report.pdf")
  file("${params.name}.params.txt")
  
  script:
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






