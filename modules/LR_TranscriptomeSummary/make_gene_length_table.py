# Module that makes a table of gencode genes and their statistical lengths
# 
#
#### READ ME ####
# The function uses a gencode fafsa file to prepare a table of genes and their statistical lengths
# 
# GenLenTab(fa_file,gen_isolen_tab = './gen_len.tsv'):
#                               inputs: 
#                                       fa_file: gencode fasta file
#                               optional inputs:
#                                       gen_len_tab: path for output table
#                               outputs:


## Import Modules ##
from pathlib import Path
import pandas as pd
import csv

## Function Stuff ##
def GenLenTab(fa_file, gen_len_stats_tab = './gen_len_stats.tsv',gen_isolen_tab='./gen_len.tsv'):
    
    # Check to see if the gene statistics file exists
    stats_file = Path(gen_len_stats_tab)
    if stats_file.is_file()==False:
        
        # If the isoform -> gene -> length file (isoinfo) from gencode does not exist, make the file
        len_file = Path(gen_isolen_tab)
        if len_file.is_file()== False:
            # Initialize Lists
            isos = []
            genes = []
            lens = []

            # Parse fafsa file for isoform, gene and length information  and append to lists
            for line in open(fa_file):
                if line.startswith('>'):
                    isos.append(line.split('|')[4].split('""')[0])
                    genes.append(line.split('|')[5].split('""')[0])
                    lens.append(line.split('|')[6].split('""')[0])
            
            # Save file as a tsv
            with open(gen_isolen_tab, 'w') as ofile:
                writer = csv.writer(ofile,delimiter='\t')
                writer.writerows(zip(isos, genes, lens))
            print('Prepared the isoform length table.')
        else:
            print('The jurkat isoform information file prepared from the gencode fafsa file exists. Skipping this step')

        # Import isolen info table #
        data = pd.read_csv(gen_isolen_tab, delimiter='\t', header=None)
        data.columns = ['isoform', 'gene', 'length']
        cut = data[['gene', 'length']]

        # Map genes to lengths and calc average, min and max. Round mean to nearest tenth. Reset indices. 
        len= cut.groupby(['gene']).length.agg(['mean', 'min', 'max'])
        len['mean'] = len['mean'].round(decimals = 1)
        len = len.reset_index(level=['gene'])

        # Change column names and save the table 
        len.columns =['gene', 'avg_len', 'min_len', 'max_len']
        len.to_csv(gen_len_stats_tab, sep="\t", index=False)
        len = pd.read_csv(gen_len_stats_tab, sep='\t')
        print('Prepared the gene length statistics table')
    else:
        len = pd.read_csv(gen_len_stats_tab, sep='\t')
        print('The gene statistics table has already been prepared. Skipping this step')

    return len

    
    