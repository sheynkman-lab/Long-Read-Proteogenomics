script
python3 prepare_reference_table_v2.py \
	-g ../../data/gencode.v35.annotation.gtf \ 
	-fa ../../data/gencode.v35.pc_transcripts.fa \ 
	--ensg_gene ../../results/PG_ReferenceTables/ensg_to_gene.tsv \
	--enst_isoname ../../results/PG_ReferenceTables/enst_to_isoname.tsv \
	--gene_ensp ../../results/PG_ReferenceTables/gene_to_ensp.tsv \
	--gene_isoname ../../results/PG_ReferenceTables/gene_to_isoname.tsv \
	--isoname_lens ../../results/PG_ReferenceTables/isoname_lens.tsv \
	--gen_lens ../../results/PG_ReferenceTables/gene_lens.tsv \
