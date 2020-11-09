makeblastdb -in a_gencode_fasta_same_prot_clustered.v35.fa -dbtype prot







blastp -db ./a_gencode_fasta_same_prot_clustered.v35.fa -query ./uniprot_reviewed_canonical_and_isoform.fasta -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send evalue bitscore length nident gaps" -out run_UN_to_GC_blast_result.tsv -num_threads 24 >> run_UN_to_GC_blast_result_stdout.txt 2>&1

