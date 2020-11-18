# Stop with info message, exit status 0
nextflow run lr_orfcalling.nf --trans_decoder false

# Stop with error message "No such file: ", the file does not exit in the location, exit status 1
nextflow run lr_orfcalling.nf --fasta https://zenodo.org/record/4278034/files/toy_for_christina.fast

# Complete successfully, exit status 0
nextflow run lr_orfcalling.nf --fasta https://zenodo.org/record/4278034/files/toy_for_christina.fasta.txt --trans_decoder true