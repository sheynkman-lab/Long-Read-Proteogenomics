#!/usr/bin/env python3
# code from making custom regions of gencode + pacbio

# %%

import subprocess
import pandas as pd
import os
import argparse
import gtfparse




def make_region_bed_for_ucsc(name, sample_gtf, reference_gtf):
    # exon bed ranges of pb
    pb = pd.read_table(sample_gtf, skiprows=1, header=None)
    pb = pb[pb[2] == 'CDS']
    # only move forward with entries with a mapped gene
    pb = pb[~pb[8].str.contains('gene_id "-";')]
    pb = pb[[0, 3, 4]]
    both = pb
    
    if reference_gtf is not None:
        gc = pd.read_table(reference_gtf, skiprows=5, header=None)
        gc = gc[gc[2] == 'exon']
        gc = gc[[0, 3, 4]]
        both = pd.concat([pb, gc])
    
    both.columns = ['chr', 'start', 'end']
    both['start'] = both['start'] - 1

    # delete after processing pb only
    # gc = pd.read_table(reference_gtf, skiprows=5, header=None)
    # gc = gc[gc[2] == 'exon']
    # gc = gc[[0, 3, 4]]

    # # concatenate pb and gc
    # both = pd.concat([pb, gc])



    fn = f'{name}_cds_ranges.bed'
    both.to_csv(fn, sep='\t', index=None, header=None)

    # sort bed
    rows = subprocess.getoutput(f'bedtools sort -i {fn}')
    fn_sorted = fn.replace('.', '_sorted.')

    # subtract 1 from start
    with open(fn_sorted, 'w') as ofile:
        ofile.write(rows)
        # for row in rows:            
        #     chr, start, end = row.split("\t")
        #     start = int(start)
        #     start = start - 1
        #     ofile.write(f"{chr}\t{start}\t{end}\n")

    # merge ranges
    # rows = subprocess.getoutput(['bedtools', 'merge', '-i', fn_sorted])
    rows = subprocess.getoutput(f'bedtools merge -i {fn_sorted}')

    # write to file
    with open(f"{name}_ucsc_multiregion.bed", "w") as ofile:
    # with open('pb_gc_ranges_for_ucsc_multregion_pb_only.bed', 'w') as ofile:
        # ofile.write('#database hg38\n#shortDesc Exons of isoforms in GENCODE v35 and Jurkat PacBio\n#padding 20\n')
        # uncomment after looking at pb only
        ofile.write(f'#database hg38\n#shortDesc Exons of isoforms in {name} PacBio CDS\n#padding 10\n')
        ofile.write(rows)

    # clean-up files
    os.remove(fn)
    os.remove(fn_sorted)



def main():
    parser = argparse.ArgumentParser("IO file locations for making region bed")
    parser.add_argument("--name", action="store", dest="name", help="name of sample - used for output file name")
    parser.add_argument("--sample_gtf", action="store", dest = "sample_gtf", help="sample gtf with cds. from make_pacbio_cds_gtf")
    parser.add_argument("--reference_gtf", action="store", dest="reference_gtf", help = "reference gtf", default="")
    results = parser.parse_args()
    make_region_bed_for_ucsc(results.name, results.sample_gtf, results.reference_gtf)


if __name__ == "__main__":
    main()


