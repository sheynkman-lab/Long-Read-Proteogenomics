bash download_pwiz_sandbox.sh
sed -i.bak 's/chambm\/pwiz-skyline-i-agree-to-the-vendor-licenses/.\/singularity_images\/pwiz_sandbox/g' \
    nextflow.config &&
    rm nextflow.config.bak
sed -i.bak 's/docker.enabled/singularity.enabled/g' nextflow.config &&
    rm nextflow.config.bak
