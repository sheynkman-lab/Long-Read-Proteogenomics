#!/usr/bin/env python3

# add rgb shading value based on the relative abundances of all pb transcripts
# of a gene

# %%

import pandas as pd
import math
import argparse

# examine all pb transcripts of a gene, determine rgb color
def calculate_rgb_shading(grp):
    """
    Examine CPM for all PB transc
    ripts of a gene and get rgb shading factor.
    """

    # rgb scaling
    rgb_scale = [
        '0,0,0', '26,0,0', '51,0,0', '77,0,0', '102,0,0',
        '128,0,0', '153,0,0', '179,0,0', '204,0,0', '230,0,0',
        '255,0,0', '255,26,26', '255,51,51', '255,77,77', '255,102,102',
        '255,128,128', '255,153,153', '255,179,179', '255,204,204', '255,230,230']
    max_cpm = grp.cpm.max()
    out_df = pd.DataFrame(columns = ['acc_full',  'pb_acc', 'cpm', 'fc', 'log2fc', 'log2fcx3', 'ceil_idx', 'rgb'])
    for i, row in grp.iterrows():
        cpm = row['cpm']
        fc = float(max_cpm) / float(cpm)
        log2fc = math.log(fc, 2) 
        log2fcx3 = log2fc * 3
        ceil_idx = math.ceil(log2fcx3)
        if ceil_idx > 19:
            ceil_idx = 19
        rgb = rgb_scale[ceil_idx] 
        out_df = out_df.append({'acc_full': row['acc_full'],
                       'pb_acc': row['pb_acc'],
                       'cpm': row['cpm'],
                       'fc': fc,
                       'log2fc': log2fc,
                       'log2fcx3': log2fcx3,
                       'ceil_idx': ceil_idx,
                       'rgb': rgb}, ignore_index=True)
    # comment out line below to return all intermediate values
    out_df = out_df[['acc_full', 'pb_acc', 'cpm', 'fc', 'rgb']]
    return out_df


def add_rgb_shading(name, bed_file):
    """
    Reads a BAM file containing CPM info to determine rgb color to use for track visualizatio

    Parameters
    ----------
    name : str 
        name of sample
    bed_file : filename
        file of bed cds to read
    """
    bed = pd.read_table(bed_file, header=None)
    bed[['gene', 'pb_acc', 'cpm']] = bed[3].str.split('|', expand=True)
    bed = bed[bed.gene != '-']
    bed = bed.rename(columns={3: 'acc_full'})
    
    # subset df to determine rgb shading
    subbed = bed[['acc_full', 'gene', 'pb_acc', 'cpm']].copy()
    subbed['cpm'] = subbed['cpm'].astype(str).astype(float)

    shaded = subbed.groupby('gene').apply(calculate_rgb_shading).reset_index()

    # include rgb into original bed12
    shaded['cpm_int'] = shaded['cpm'].apply(lambda x: str(round(x)).split('.')[0])
    shaded['new_acc_full'] = shaded['gene'] + '|' + shaded['pb_acc'] + '|' + shaded['cpm_int'].astype(str)

    # join in the rgb data and new accession
    bed_shaded = pd.merge(bed, shaded, how='left', on='acc_full')
    bed_shaded = bed_shaded[[0, 1, 2, 'new_acc_full', 4, 5, 6, 7, 'rgb', 9, 10, 11]]

    with open(f'{name}_cds_shaded.bed12', 'w') as ofile:
        ofile.write('track name=pacbio_cds_w_rgb_shade itemRgb=On\n')
        bed_shaded.to_csv(ofile, sep='\t', index=None, header=None)

    
    


def main():
    parser = argparse.ArgumentParser("IO file locations for making region bed")
    parser.add_argument("--name", action="store", dest="name", help="name of sample - used for output file name")
    parser.add_argument("--bed_file", action="store", dest = "bed_file", help="sample bed with cds")
    results = parser.parse_args()
    add_rgb_shading(results.name, results.bed_file)




if __name__ == "__main__":
    main()