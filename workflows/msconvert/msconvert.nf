/*--------------------------------------------------
Peptide Track Visualization
 * Makes peptide track for UCSC Genome Browser 
---------------------------------------------------*/

if (!params.mass_spec == false) {
   if (!params.mass_spec.endsWith("tar.gz")) {
      ch_mass_spec_raw = Channel.fromPath("${params.mass_spec}/*.raw")
   } else {
      if (params.mass_spec.endsWith("tar.gz")){
         ch_mass_spec_raw_mzml_tar_gz = Channel.value(file(params.mass_spec))
      }
   }
} else {
   ch_mass_spec_raw = Channel.from("no mass spec")
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
         script:
            """
            tar xvzf $raw_mzml_tar_gz
            """
      }
   }
}

process mass_spec_raw_convert{
   tag "${raw_file}"
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