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
    out_df = out_df[['acc_full','rgb']]
    return out_df


def add_rgb_shading_cpm(name, bed,split_size):
    """
    Reads a BAM file containing CPM info to determine rgb color to use for track visualizatio

    Parameters
    ----------
    name : str 
        name of sample
    bed_file : filename
        file of bed cds to read
    """
    
    
    # subset df to determine rgb shading
    if split_size==3:
        subbed = bed[['acc_full', 'gene', 'pb_acc', 'cpm']].copy()
    elif split_size==4:
        subbed = bed[['acc_full', 'gene', 'pb_acc','pclass', 'cpm']].copy()
    subbed['cpm'] = subbed['cpm'].astype(str).astype(int)

    shaded = subbed.groupby('gene').apply(calculate_rgb_shading).reset_index(drop=True)

    # include rgb into original bed12
    bed_shaded = pd.merge(bed, shaded, how='left', on='acc_full')
    
    # shaded['cpm_int'] = shaded['cpm'].apply(lambda x: str(round(x)).split('.')[0])
    if split_size==3:
        bed_shaded['new_acc_full'] = bed_shaded.apply(lambda row: '|'.join([row['gene'],row['pb_acc'],str(row['cpm'])]),axis=1)
    if split_size==4:
        bed_shaded['new_acc_full'] = bed_shaded.apply(lambda row: '|'.join([row['gene'],row['pb_acc'],row['pclass'],str(row['cpm'])]),axis=1)


    # join in the rgb data and new accession
    
    bed_shaded = bed_shaded[['chrom', 'chromStart', 'chromStop', 'new_acc_full', 'score', 'strand', 'thickStart', 'thickEnd', 'rgb', 'blockCount', 'blockSizes', 'blockStarts']]

    with open(f'{name}_cds_shaded_cpm.bed12', 'w') as ofile:
        ofile.write(f'track name={name}_cds_w_cpm_rgb_shade itemRgb=On\n')
        bed_shaded.to_csv(ofile, sep='\t', index=None, header=None)


def add_rgb_shading_pclass(name,bed):
    pclass_shading_dict = {
        'pFSM':'100,165,200',
        'pNIC':'111,189,113',
        'pNNC':'232,98,76',
        'pISM':'248,132,85'
    }
    bed['rgb'] = bed['pclass'].map(pclass_shading_dict).fillna('0,0,0')
    filter_names = ['chrom','chromStart','chromStop','acc_full','score','strand','thickStart','thickEnd','rgb','blockCount','blockSizes','blockStarts']
    bed = bed[filter_names]
    with open(f'{name}_cds_shaded_protein_class.bed12', 'w') as ofile:
        ofile.write(f'track name={name}_cds_w_pclass_rgb_shade itemRgb=On\n')
        bed.to_csv(ofile, sep='\t', index=None, header=None)
    
def add_rgb_shading(name, bed_file):
    bed_names = ['chrom','chromStart','chromStop','acc_full','score','strand','thickStart','thickEnd','itemRGB','blockCount','blockSizes','blockStarts']
    bed = pd.read_table(bed_file, names=bed_names)
    split_size=len(bed.loc[0,'acc_full'].split('|'))
    if split_size==3:
        bed[['gene', 'pb_acc', 'cpm']] = bed['acc_full'].str.split('|', expand=True)
    if split_size==4:
        bed[['gene', 'pb_acc','pclass', 'cpm']] = bed['acc_full'].str.split('|', expand=True)
    bed = bed[bed.gene != '-']

    add_rgb_shading_cpm(name, bed.copy(), split_size)
    if split_size==4:
        add_rgb_shading_pclass(name,bed)


def main():
    parser = argparse.ArgumentParser("IO file locations for making region bed")
    parser.add_argument("--name", action="store", dest="name", help="name of sample - used for output file name")
    parser.add_argument("--bed_file", action="store", dest = "bed_file", help="sample bed with cds")
    results = parser.parse_args()
    add_rgb_shading(results.name, results.bed_file)




if __name__ == "__main__":
    main()
# %%
