// Channel
//     .value(file(params.toml))
//     .ifEmpty { error "Cannot find toml file for parameter --toml: ${params.toml}" }
//     .set { ch_toml } 

// Channel
//     .value(file(params.orf_calls))
//     .ifEmpty { error "Cannot find gtf file for parameter --gencode_gtf: ${params.orf_calls}" }
//     .set { ch_orf_calls } 


Channel
    .value(file(params.orf_fasta))
    .ifEmpty { error "Cannot find gtf file for parameter --orf_fasta: ${params.orf_fasta}" }
    .set { ch_orf_fasta }




Channel
    .fromPath("${params.mass_spec}/*.raw")
    .set{ch_mass_spec_raw}

Channel
    .fromPath("${params.mass_spec}/*.{mzml,mzML}")
    .set{ch_mass_spec_mzml}


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

ch_mass_spec_combined = ch_mass_spec_mzml.concat(ch_mass_spec_converted)

process metamorpheus_with_sample_specific_database{
    tag " $mass_spec"
    publishDir "${params.outdir}/metamorpheus/", mode: 'copy'

    input:
        // file(orf_calls) from ch_orf_calls
        file(orf_fasta) from ch_orf_fasta
        // file(toml) from ch_toml
        file(mass_spec) from ch_mass_spec_combined.collect()

    output:
        file("toml/*")
        file("search_results/*")
        file("search_reseults/Task1SearchTask/AllPeptides.psmtsv") into ch_sample_specific_peptides
    
    script:
        """
        dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings settings 
        dotnet /metamorpheus/CMD.dll -d $orf_fasta -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results
        """
}

// process metamorpheus_with_gencode_database{
//     tag " $mass_spec"
//     publishDir "${params.outdir}/metamorpheus/", mode: 'copy'

//     input:
//         // file(orf_calls) from ch_orf_calls
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



