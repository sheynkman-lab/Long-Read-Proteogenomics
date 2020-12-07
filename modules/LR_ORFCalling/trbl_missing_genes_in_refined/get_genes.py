# %%
from Bio import SeqIO
from collections import defaultdict

gene_pb = defaultdict(list)

for rec in SeqIO.parse('jurkat_orf_refined.fasta', 'fasta'):
	gene = rec.description.split('=')[1]
	pb = rec.id
	gene_pb[gene].append(pb)

# %%
	
	
with open('genes_in_refined.tsv', 'w') as ofile:
	for gene, pbs in gene_pb.items():
		ofile.write(gene + '\t' + ','.join(pbs) + '\n')