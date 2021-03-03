nextflow sqanti3.nf \
--name test_sqanit3 \
--gencode_fasta /mnt/shared/ubuntu/session_data/data/test_data/gencode.v35.pc_transcripts.chr22.fa \
--gencode_gtf   /mnt/shared/ubuntu/session_data/data/test_data/gencode.v35.annotation.chr22.gtf \
--sample_gtf    /mnt/shared/ubuntu/session_data/data/test_data/jurkat.collapsed.chr22.gff \
--sample_fasta  /mnt/shared/ubuntu/session_data/data/test_data/jurkat.collapsed.chr22.fasta \
--fl_count      /mnt/shared/ubuntu/session_data/data/jurkat.collapsed.abundance.txt 