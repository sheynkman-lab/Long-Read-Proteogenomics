nextflow orf_calling.nf \
--name orf-call \
--orf_coord /mnt/shared/ubuntu/session_data/data/test_data/jurkat_cpat.ORF_prob.chr22.tsv \
--gencode_gtf /mnt/shared/ubuntu/session_data/data/test_data/gencode.v35.annotation.chr22.gtf \
--sample_gtf /mnt/shared/ubuntu/session_data/data/test_data/jurkat.collapsed.chr22.2.gff \
--pb_gene /mnt/shared/ubuntu/session_data/data/test_data/pb_to_gene.tsv \
--classification /mnt/shared/ubuntu/session_data/data/test_data/jurkat_classification.txt \
--sample_fasta /mnt/shared/ubuntu/session_data/data/test_data/jurkat.collapsed.chr22.fasta