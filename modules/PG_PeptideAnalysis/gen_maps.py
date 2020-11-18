# Module that extracts tables of 1. ENSG-genename
#                                2. ENST-isoname
#                                3. genename-protein

#### READ ME ####
# The function GenMap maps ENSGs to genenames, ENSTs to transcript names
# and gene names to ENSPs
# 
# GenMap(gtf_file, create_files):
#                               inputs: 
#                                       gtf file: gtf file from gencode
#                               optional inputs:
#                                       create_files = all or a combination of ensg_gene
#                                                      enst_isoname and gene_ensp 
#                                                      
#                                                      all is the default values
#                                                      ensg_to_gene: makes a tsv file of ensg -> genename
#                                                      enst_to_trans: makes a tsv file of enst -> transcript name
#                                                      ensp_to_gene: makes a tsv file of gene -> ensp 
#                                                      trans_to_gene: makes a tsv tile of trans -> gene 
#
#                                                      usage:
#                                                          Generate all files: GenMap(gtf_file)
#                                                          Make ensg_gene and gene_ensp files: 
#                                                           GenMap(gtf_gile, ["ensg_to_gene", "ensp_to_gene"])  
#  

## Import Modules ##
from collections import defaultdict

def GenMap(gtf_file, create_files="all"):
    genes = {} # ENSG -> <genename>
    isos = {} # ENST -> <iso-name>
    ensps = defaultdict(set) # gene_name -> <set(ENSPs)>
    trans = defaultdict(set) # transcript name -> gene_name 

    for line in open(gtf_file):
        if line.startswith('#'):
            pass
        else:
            wds = line.split('\t')
            cat = wds[2]
            if cat in ['transcript']:
                ensg = line.split('gene_id "')[1].split('"')[0].replace('-', '_')
                gene = line.split('gene_name "')[1].split('"')[0].replace('-', '_')
                enst = line.split('transcript_id "')[1].split('"')[0].replace('-', '_')
                transcript = line.split('transcript_name "')[1].split('"')[0].replace('-', '_')
                if "ensg_to_gene" in create_files or create_files == "all":
                    genes[ensg] = gene
                if "enst_to_trans" in create_files or create_files == "all":
                    isos[enst] = transcript
                if "trans_to_gene" in create_files or create_files == 'all':
                    trans[gene].add(transcript)
                if "ensp_to_gene" in create_files or create_files == "all":
                    if 'transcript_type "protein_coding' in line:
                        gen = line.split('gene_name "')[1].split('"')[0].replace('-', '_')
                        ensp = line.split('protein_id "')[1].split('"')[0].replace('-', '_')
                        ensps[gen].add(ensp)
    
    if "ensg_to_gene" in create_files or create_files == "all":
        with open('ensg_to_gene.tsv', 'w') as ofile:
            for ensg, gene in genes.items():
                ofile.write(ensg + '\t' + gene + '\n')
        print("File ensg_to_gene.tsv prepared")

    if "enst_to_trans" in create_files or create_files == "all":
        with open('enst_to_trans.tsv', 'w') as ofile:
            for enst, isoname in isos.items():
                ofile.write(enst + '\t' + isoname + '\n')
        print("File enst_to_trans.tsv prepared")
    
    if "ensp_to_gene" in create_files or create_files == "all":
        with open('ensp_to_gene.tsv', 'w') as ofile:
            for gen, ensp_set in ensps.items():
                for ensp in ensp_set:
                    ofile.write(gen + '\t' + ensp + '\n')
        print("File ensp_to_gene.tsv prepared")
    
    if "trans_to_gene" in create_files or create_files == "all":
        with open('trans_to_gene.tsv', 'w') as ofile:
            for gene, transcript_set in trans.items():
                for transcript in transcript_set:
                    ofile.write(gene + '\t' + transcript + '\n')
        print("File trans_to_gene.tsv prepared")
    
    
    
    

    

