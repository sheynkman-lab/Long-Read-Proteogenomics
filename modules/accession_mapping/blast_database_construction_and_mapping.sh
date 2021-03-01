# inter-database isoform mapping with blast

# install blast
conda install -c bioconda blast

# note - ensure that the databases going in have unique protein sequences

# make blast databases
makeblastdb -in gencode.fasta -dbtype prot
makeblastdb -in uniprot.fasta -dbtype prot
makeblastdb -in pacbio.fasta -dbtype prot

# do blast searches
# pairwise between databases (un to gc, gc to un, pb to gc, gc to pb, pb to un, un to pb)
blastp -db ${blast_db} -query ${query_fasta} -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send evalue bitscore length nident gaps" -out query_to_database_blast_result.tsv -num_threads 40 >> query_to_database_stdout.txt 2>&1

# find "at-length" isoform matches between the databases
find_at_length_blast_matches.py




