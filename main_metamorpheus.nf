// Channel
//     .value(file(params.toml))
//     .ifEmpty { error "Cannot find toml file for parameter --toml: ${params.toml}" }
//     .set { ch_toml } 

// Channel
//     .value(file(params.orf_calls))
//     .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.orf_calls}" }
//     .set { ch_orf_calls } 


// TODO - is orf_fasta the refined database? If so, can we call this parameter --refined_protein_db or somethign like that?
// TODO cont. - orf typically means the nt sequence underlying the protein



if (!params.sqanti_fasta) exit 1, "Cannot find any sqanti_fasta file for parameter --sqanti_fasta: ${params.sqanti_fasta}"
ch_orf_fasta = Channel.value(file(params.sqanti_fasta))

if (!params.orf_meta) exit 1, "Cannot find any sqanti_fasta file for parameter --orf_meta: ${params.orf_meta}"
ch_orf_meta = Channel.value(file(params.orf_meta))

Channel
    .fromPath("${params.mass_spec}/*.raw")
    .set{ch_mass_spec_raw}

Channel
    .fromPath("${params.mass_spec}/*.{mzml,mzML}")
    .set{ch_mass_spec_mzml}


if (!params.rescue_resolve_toml) exit 1, "Cannot find file for parameter --rescue_resolve_toml: ${params.rescue_resolve_toml}"
ch_rr_toml = Channel.value(file(params.rescue_resolve_toml))

ch_metamorpheus_toml = Channel.value(file(params.metamorpheus_toml))


process mass_spec_raw_convert{
    publishDir "${params.outdir}/raw_convert/", mode: 'copy'

    input:
        file(raw_file) from ch_mass_spec_raw

    output:
        file("*") into ch_mass_spec_converted

    script:
        """
        wine msconvert $raw_file --filter "peakPicking true 1-"
        """
}

// TODO - can you explain to me why this was needed?
ch_mass_spec_combined = ch_mass_spec_mzml.concat(ch_mass_spec_converted)
ch_mass_spec_combined.into{
    ch_mass_spec_for_pacbio_rescue_resolve
    ch_mass_spec_for_pacbio_normal
}

process metamorpheus_with_sample_specific_database{
    tag " $mass_spec"
    label 'metamorpheus'
    publishDir "${params.outdir}/metamorpheus/", mode: 'copy'

    input:
        // file(orf_calls) from ch_orf_calls
        file(orf_fasta) from ch_orf_fasta
        file(toml_file) from ch_metamorpheus_toml
        file(mass_spec) from ch_mass_spec_for_pacbio_normal.collect()

    output:
        file("*")
    
    script:
        def toml = toml_file.name != 'NO_TOML_FILE' ? "$toml_file" : 'toml/SearchTask.toml'
        """
        echo $toml > toml_given.txt
        echo $toml_file > toml_used.txt
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings settings 
        dotnet /metamorpheus/CMD.dll -d $orf_fasta -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results
        """
}

// process metamorpheus_with_sample_specific_database_rescue_resolve{
//     tag " $mass_spec"
//     publishDir "${params.outdir}/metamorpheus/rescue_resolve", mode: 'copy'

//     input:
//         // file(orf_calls) from ch_orf_calls
//         file(orf_fasta) from ch_orf_fasta
//         // file(toml) from ch_toml
//         file(mass_spec) from ch_mass_spec_for_pacbio_rescue_resolve.collect()
//         file(toml) from ch_rr_toml
//         file(orf_meta) from ch_orf_meta

//     output:
//         file("toml/*")
//         file("search_results/Task1SearchTask/All*")
//         file("search_results/Task1SearchTask/prose.txt")
//         file("search_results/Task1SearchTask/results.txt")
//         file("search_results/Task1SearchTask/AllPeptides.${params.name}.rescue_resolve.psmtsv")
//         file("search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.rescue_resolve.tsv")
    
//     script:
//         """
//         dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings settings 
//         dotnet /metamorpheus/CMD.dll -d $orf_fasta -s $mass_spec -t $toml -v normal --mmsettings settings -o ./search_results --orf $orf_meta --cpm 25
//         mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.${params.name}.rescue_resolve.psmtsv
//         mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.${params.name}.rescue_resolve.tsv
//         """
// }

// process metamorpheus_with_gencode_database{
//     tag " $mass_spec"
//     publishDir "${params.outdir}/metamorpheus/", mode: 'copy'

//     input:
//         // file(orf_calls) from ch_orf_calls
           // TODO - is the below variable specific enough? e.g., ch_gencode_translations_fasta
//         file(gencode_fasta) from ch_gencode_fasta
//         // file(toml) from ch_toml
//         file(mass_spec) from ch_mass_spec_combined.collect()

//     output:
//         file("*")
//         file("AllPeptides.psmtsv") into ch_gencode_peptides
    
//     script:
//         """
//         dotnet /metamorpheus/CMD.dll -g -o ./ --mmsettings settings
//         dotnet /metamorpheus/CMD.dll -d $gencode_fasta -s $mass_spec -t SearchTask.toml -v normal --mmsettings settings
//         """
// }



