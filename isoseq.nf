
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

  bamtools filter -tag 'rq':'>=0.90' -in $sample_css -out filtered.$sample_css 
  # create an index
  pbindex filtered.$sample_ccs

  
  lima --isoseq --dump-clips --peek-guess -j ${task.cpus} filtered.$sample_ccs $primers_fasta ${params.name}.demult.bam
  isoseq3 refine --require-polya ${params.name}.demult.NEB_5p--NEB_3p.bam $primers_fasta ${params.name}.flnc.bam

  # clustering of reads, can only make faster by putting more cores on machine (cannot parallelize)
  isoseq3 cluster ${params.name}.flnc.bam ${params.name}.clustered.bam --verbose --use-qvs

  # align reads to the genome, takes few minutes (40 core machine)
  pbmm2 align $gencode_fasta ${params.name}.clustered.hq.bam ${params.name}.aligned.bam --preset ISOSEQ --sort -j ${task.cpus} --log-level INFO

  # collapse redundant reads
  isoseq3 collapse ${params.name}.aligned.bam ${params.name}.collapsed.gff
  """
}