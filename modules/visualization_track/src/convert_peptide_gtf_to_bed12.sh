bname="test_peptides"
gtfToGenePred ${bname}.gtf ${bname}.genePred
genePredToBed ${bname}.genePred ${bname}.bed12

# add rgb to colorize specific peptides 
echo 'track name=peptide_w_specificity itemRgb=On' > ${bname}_rgb.bed12
cat ${bname}.bed12 | awk '{ if ($4 ~ /-1$/) {$9="0,102,0"; print $0} else {$9="0,51,0"; print $0} }' >> ${bname}_rgb.bed12