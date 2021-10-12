/*--------------------------------------------------
Peptide Track Visualization
 * Makes peptide track for UCSC Genome Browser 
---------------------------------------------------*/

if (!params.sample_gtf) exit 1, "Cannot find any file for parameter --gtf: ${params.sample_gtf}"
ch_sample_gtf = Channel.value(file(params.sample_gtf))

if (!params.reference_gtf) exit 1, "Cannot find any file for parameter --reference_gtf: ${params.reference_gtf}"
ch_reference_gtf = Channel.value(file(params.reference_gtf))

if (!params.peptides) exit 1, "Cannot find any file for parameter --peptides: ${params.peptides}"
ch_peptides = Channel.value(file(params.peptides))

if (!params.sample_fasta) exit 1, "Cannot find any file for parameter --sample_fasta: ${params.sample_fasta}"
ch_sample_fasta = Channel.value(file(params.sample_fasta))

if (!params.pb_gene) exit 1, "Cannot find any file for parameter --pb_gene: ${params.pb_gene}"
ch_pb_gene = Channel.value(file(params.pb_gene))

if (!params.peptides) exit 1, "Cannot find any file for parameter --peptides: ${params.gene_isoname}"
ch_gene_isoname = Channel.value(file(params.gene_isoname))


process make_peptide_gtf_file{
    publishDir "${params.outdir}"
    input:
        file(sample_gtf) from ch_sample_gtf
        file(reference_gtf) from ch_reference_gtf
        file(peptides) from ch_peptides
        file(sample_fasta) from ch_sample_fasta
        file(pb_gene) from ch_pb_gene
        file(gene_isoname) from ch_gene_isoname
    
    output:
        file("${params.name}_peptides.gtf") into ch_peptide_gtf
    script:
        """
        make_peptide_gtf_file.py \
        --name ${params.name} \
        --sample_gtf $sample_gtf \
        --reference_gtf $reference_gtf \
        --peptides $peptides \
        --pb_gene $pb_gene \
        --gene_isoname $gene_isoname \
        --refined_fasta $sample_fasta
        """
}
process generate_peptide_bed_file{
    publishDir "${params.outdir}"
    input:
        file(peptide_gtf) from ch_peptide_gtf
    output:
        file("*.bed12")
    script:
        """
        gtfToGenePred ${params.name}_peptides.gtf ${params.name}_peptides.genePred
        genePredToBed ${params.name}_peptides.genePred ${params.name}_peptides.bed12
        # add rgb to colorize specific peptides 
        finalize_peptide_bed.py \
        --bed ${params.name}_peptides.bed12 \
        --name ${params.name}
        """
        

}