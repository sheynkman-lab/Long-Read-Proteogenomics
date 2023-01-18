sed -i.bak 's/.\/singularity_images\/pwiz_sandbox/chambm\/pwiz-skyline-i-agree-to-the-vendor-licenses/g' \
    nextflow.config &&
    rm nextflow.config.bak
sed -i.bak 's/singularity.enabled/docker.enabled/g' nextflow.config &&
    rm nextflow.config.bak
