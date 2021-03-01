# commands gloria previously ran to process pb transcripts
# through sqanti
# should work here too

# input:
# polished transcripts from iso-seq (jurkat.collapsed.fasta)
# gencode annotation (gencode.v35.annotation.gtf)
# genome, only 'canonical' chromosomes  (hg38_canon.fa)
# isoform counts (jurkat.collapsed.abundance.txt)

# output of interest:
# isoform classification table

python sqanti3_qc.py jurkat.collapsed.fasta gencode.v35.annotation.gtf hg38_canon.fa -o jurkat -d SQANTI3_out_v2/ --fl_count jurkat.collapsed.abundance.txt -n8

# note - sqanti3.py does minimap alignment of pacbio transcripts against hg38

