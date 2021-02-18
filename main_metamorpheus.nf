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
    .fromPath("${params.mass_spec}/*.{raw,mzml,mzML}")
    .set{ch_mass_spec}

for(file in ch_mass_spec){
    if(file.endsWith(".raw")){
        println(file)
    }
}
//     
// Channel
//     .fromPath('data/test_metamorpheus',type:'dir',relative:true)
//     .view()

// print("\n*********************\n")
// print(ch_mass_spec)

// process raw_convert{
//     publishDir "${params.outdir}/raw_convert/", mode: 'copy'
//     input:
//         file(raw_file) from ch_mass_spec
//     output:
//         file("*") into ch_mass_spec_mzml
//     script:
//         if(raw_file.endsWith(".raw")){
//             """
//             wine msconvert $raw_file --filter "peakPicking true 1-"
//             """
//         }
//         else{
//             """
//             cp $raw_file ${params.name}.$raw_file
//             """
//         }
    

}

// process metamorpheus{
//     tag "$orf_calls, $orf_fasta, $mass_spec_fraction"
//     publishDir "${params.outdir}/metamorpheus/", mode: 'copy'

//     input:
//         // file(orf_calls) from ch_orf_calls
//         file(orf_fasta) from ch_orf_fasta
//         // file(toml) from ch_toml
//         file(mass_spec_fraction) from ch_mass_spec_mzml

//     output:
//         file("*")
    
//     script:


//         """
//         dotnet /metamorpheus/CMD.dll -g -o ./ --mmsettings settings
//         sed -i 's/false/true/g' settings/settings.toml

//         dotnet /metamorpheus/CMD.dll -d $orf_fasta -s $mass_spec_fraction -t SearchTask.toml -v normal --mmsettings settings
//         """
//         // """
//         // dotnet /metamorpheus/CMD.dll --test -v minimal -o metamorpheus
//         // """

//         // """
//         // echo "y" | dotnet /metamorpheus/CMD.dll \
//         // -d $orf_fasta \
//         // -s $mass_spec_fraction \
//         // -t $toml \
//         // -v normal \
//         // """
// } 
