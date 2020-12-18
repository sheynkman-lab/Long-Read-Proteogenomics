
python bin/refine_orf.py \
--orfs                 /mnt/shared/ubuntu/session_data/data/results/orf_calling/jurkat.collapsed.ORF_called.tsv \
--pb_fasta             /mnt/shared/ubuntu/session_data/data/jurkat.collapsed.fasta \
--redundant            /mnt/shared/ubuntu/session_data/data/results/refine/no_filter/jurkat.collapsed_redundant_accessions.txt \
--combined_tsv         /mnt/shared/ubuntu/session_data/data/results/refine/no_filter/jurkat.collapsed_orf_combined.tsv \
--combined_fasta       /mnt/shared/ubuntu/session_data/data/results/refine/no_filter/jurkat.collapsed_orf_combined.fasta \
--agg_tsv              /mnt/shared/ubuntu/session_data/data/results/refine/no_filter/jurkat.collapsed_orf_aggregated.tsv \
--agg_fasta            /mnt/shared/ubuntu/session_data/data/results/refine/no_filter/jurkat.collapsed_orf_aggregated.fasta \
--protein_coding_only  no \
--protein_coding_genes /mnt/shared/ubuntu/session_data/data/results/protein_coding_genes.txt \
--coding_score_cutoff 0.0


