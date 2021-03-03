docker \
   run \
   -it \
   --mount target=/data,source=/mnt/shared/ubuntu/session_data/data,type=bind \
   milescsmith/sqanti3 \
   sqanti3_qc \
        data/test_data/jurkat.collapsed.chr22.fasta \
        data/test_data/gencode.v35.annotation.chr22.gtf \
        data/test_data/gencode.v35.pc_transcripts.chr22.fa \
        -o jurkat \
        -d SQANTI3_out/ \
        --fl_count data/jurkat.collapsed.abundance.txt \