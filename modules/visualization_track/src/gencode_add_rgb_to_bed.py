#!/usr/bin/env python3

import pandas as pd 
import argparse

def add_rgb_colors(bed_file, rgb, version):
    bed_names = ['chrom','chromStart','chromStop','acc_full','score','strand','thickStart','thickEnd','itemRGB','blockCount','blockSizes','blockStarts']
    bed = pd.read_table(bed_file, names=bed_names)
    bed['rgb'] = rgb
    filter_names = ['chrom','chromStart','chromStop','acc_full','score','strand','thickStart','thickEnd','rgb','blockCount','blockSizes','blockStarts']
    bed = bed[filter_names]
    with open('gencode_shaded.bed12', 'w') as ofile:
        ofile.write(f'track name="GENCODE version {version}" itemRgb=On\n')
        bed.to_csv(ofile, sep='\t', index=None, header=None)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gencode_bed", action="store", dest = "gencode_bed",)
    parser.add_argument("--rgb", action="store", dest="rgb", default="0,0,140")
    parser.add_argument("--version",action="store",dest="version", default="V35")
    args = parser.parse_args()
    add_rgb_colors(args.gencode_bed, args.rgb, args.version)

if __name__ == "__main__":
    main()