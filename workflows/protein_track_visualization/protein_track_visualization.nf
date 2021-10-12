if (!params.sample_gtf) exit 1, "Cannot find any gtf file for parameter --gtf: ${params.sample_gtf}"
ch_sample_gtf = Channel.value(file(params.sample_gtf))

if (!params.reference_gtf) exit 1, "Cannot find any gtf file for parameter --reference_gtf: ${params.reference_gtf}"
ch_reference_gtf = Channel.value(file(params.reference_gtf))

process convert_gtf_to_bed{
    publishDir "${params.outdir}"
    input:
        file(gtf) from ch_sample_gtf
    output:
        file("${params.name}.bed12") into ch_bed
    script:
    """
    gtfToGenePred $gtf ${params.name}.genePred
    genePredToBed ${params.name}.genePred ${params.name}.bed12
    """
}
process add_colors_to_track{
    publishDir "${params.outdir}"
    input:
        file(bed) from ch_bed
    output:
        file("*")
    script:
    """
    track_add_rgb_colors_to_bed.py \
    --name ${params.name} \
    --bed_file $bed
    """
}

process multiregion{
    publishDir "${params.outdir}"
    input:
        file(sample_gtf) from ch_sample_gtf
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